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


"""
if __name__ == "__main__":
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
