"""Module for processing image files."""
import matplotlib.image as mpimg
from scipy import fftpack, ndimage
import numpy as np
from matplotlib import pyplot as plt

# from scipy.ndimage import rotate
import scipy


class Image(object):
    """Module for processing image files."""

    def __init__(self, values=[], pbc_params=[]):
        """Initialize the module."""
        # pbc_params for 2D image: a,b and gamma
        # 3D : a,b,c, alpha, beta, gamma

        self.values = values  # XY numpy array
        self.pbc_params = pbc_params

    def zoom_interp_2d(self, FFT_image, zoom_factor=40, interpol_factor=1):
        """Zoom and interpolate an image."""
        zoom_size = (FFT_image.shape[0] / zoom_factor) / 2
        window_size_x = FFT_image.shape[0]
        window_size_y = FFT_image.shape[1]

        if np.mod(FFT_image.shape[0] / zoom_factor, 2) == 0:
            F2_zoomed = FFT_image[
                int(window_size_x / 2 - zoom_size) : int(
                    window_size_x / 2 + zoom_size
                ),
                int(window_size_y / 2 - zoom_size) : int(
                    window_size_y / 2 + zoom_size
                ),
            ]
        else:
            F2_zoomed = FFT_image[
                int(window_size_x / 2 - zoom_size) : int(
                    window_size_x / 2 + 1 + zoom_size
                ),
                int(window_size_y / 2 - zoom_size) : int(
                    window_size_y / 2 + 1 + zoom_size
                ),
            ]

        return ndimage.zoom(F2_zoomed, interpol_factor)

    def fourier_transform2D(
        self,
        zoom_factor=30,
        interpol_factor=1,
        use_crop=True,
        pad_bright_spot=True,
    ):
        """Make 2D FT."""
        if use_crop:

            g1 = fftpack.fft2(
                (self.crop_square().rgb_to_gray().values[:, :, 0])
            )
        else:
            g1 = fftpack.fft2((self.rgb_to_gray().values[:, :, 0]))
        g2 = np.abs((fftpack.fftshift((g1))))
        g2 = np.abs((fftpack.fftshift((g1))))
        g3 = self.zoom_interp_2d(
            g2, zoom_factor=zoom_factor, interpol_factor=interpol_factor
        )
        if pad_bright_spot:

            g3 = g3 / np.max(g3)
            x = int(g3.shape[0] / 2)
            y = int(g3.shape[1] / 2)
            g3[x, y] = 0.15

        print("g3 shape", x, y)
        # Zoom in to see the FFT more clearly.

        return Image(values=g3)
        # ft_arr=
        # return ft_arr

    @staticmethod
    def from_file(path):
        """Load image frim file."""
        im = mpimg.imread(path)
        return Image(values=im)
        pass

    def rgb_to_gray(self):
        """Make RGB to Grey scale image."""
        img = np.array(self.values)
        # grayImage = np.zeros(img.shape)
        R = np.array(img[:, :, 0])
        G = np.array(img[:, :, 1])
        B = np.array(img[:, :, 2])

        R = R * 0.299
        G = G * 0.587
        B = B * 0.114

        Avg = R + G + B

        for i in range(3):
            img[:, :, i] = Avg

        return Image(values=img)

    def crop_square(self, size=None):
        """Crop an image."""
        if size is None:
            size = min(self.values.shape[0], self.values.shape[1])
        img_cropped = self.values.copy()[:size, :size, :]

        return Image(values=img_cropped)

    def save(self, filename):
        """Save an image."""
        # if size is None:
        # from matplotlib import pyplot as plt

        plt.imshow(self.values, interpolation="nearest")
        plt.savefig(filename)
        plt.close()

    def black_and_white(self, threshold=127):
        """Make black and white image."""
        bw = np.asarray(self.values).copy()
        bw[bw < threshold] = 0
        bw[bw >= threshold] = 255
        return Image(values=bw)

    def gaussian_filter(self, sigma=10):
        """Apply Gaussian filter to an image."""
        sigmax = sigma
        sigmay = sigma

        cy, cx = int(self.values.shape[0] / 2), int(self.values.shape[1] / 2)
        x = np.linspace(0, self.values.shape[0], self.values.shape[0])
        y = np.linspace(0, self.values.shape[1], self.values.shape[1])
        X, Y = np.meshgrid(x, y)
        gmask = np.exp(-(((X - cx) / sigmax) ** 2 + ((Y - cy) / sigmay) ** 2))
        new_values = self.values * gmask
        return Image(values=new_values)

    def rotate(self, angle=45):
        """Rotate an image."""
        rot = scipy.ndimage.rotate(self.values, angle)
        return Image(values=rot)

    @staticmethod
    def get_blob_angles(
        filename="image.png",
        tol=0.005,
        tmp_path="tmp.png",
        param1=100,
        param2=255,
        s1=3,
        s2=20,
        rad=10,
        color="r",
        out_name="final.png",
    ):
        """Get angles between atoms."""
        try:
            import cv2
        except ImportError as error:
            print(error.__class__.__name__ + ": " + error.message)
        # https://www.geeksforgeeks.org/
        # white-and-black-dot-detection-using-opencv-python/
        plt.clf()
        im = Image.from_file(filename)
        vals = im.black_and_white(tol).values
        plt.imshow(vals, cmap="jet", interpolation="none")
        plt.axis("off")
        plt.savefig(tmp_path)
        plt.close()
        path = tmp_path
        # reading the image in grayscale mode
        gray = cv2.imread(path, 0)
        # reading the image in grayscale mode
        th, threshed = cv2.threshold(
            gray, param1, param2, cv2.THRESH_BINARY | cv2.THRESH_OTSU
        )
        # findcontours
        cnts = cv2.findContours(
            threshed, cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE
        )[-2]

        # filter by area
        xcnts = []

        for cnt in cnts:
            if s1 < cv2.contourArea(cnt) < s2:
                xcnts.append(cnt)
        canvas = gray
        cv2.drawContours(canvas, xcnts, -1, (0, 255, 0), 1)
        fig, ax = plt.subplots()
        plt.imshow(canvas)
        plt.savefig(out_name)
        plt.close()
        my_coords = []
        for i in xcnts:

            circle1 = plt.Circle((i[0][0][0], i[0][0][1]), rad, color=color)
            my_coords.append([i[0][0][0], i[0][0][1]])
            ax.add_patch(circle1)
        #

        # printing output
        print("\nDots number: {}".format(len(xcnts)))

        def angle_between(p1, p2):
            ang1 = np.arctan2(*p1[::-1])
            ang2 = np.arctan2(*p2[::-1])
            return np.rad2deg((ang1 - ang2) % (2 * np.pi))

        all_angle_data = []
        for i in my_coords:
            for j in my_coords:
                if i != j:
                    all_angle_data.append([i, j, angle_between(i, j)])
                    print(i, j, angle_between(i, j))
        return all_angle_data


"""
if __name__ == "__main__":
    from jarvis.db.figshare import make_stm_from_prev_parchg
    make_stm_from_prev_parchg()
    Image.get_blob_angles(filename='stm_image.png')
    p = "JVASP-667_neg.jpg"
    im = Image.from_file(p)

    plt.imshow(
        im.fourier_transform2D(use_crop=True, zoom_factor=50)
        .rotate(angle=0)
        .black_and_white(threshold=0.05)
        .values,
        cmap="Greys",
    )
"""
