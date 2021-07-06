"""Module to simulate STM with Tershoff-Hamann approach."""
# Reference: https://www.nature.com/articles/s41597-021-00824-y
import matplotlib.pyplot as plt
from jarvis.io.vasp.outputs import Chgcar
import numpy as np
import matplotlib.transforms as mtransforms
import scipy
from jarvis.core.image import Image


class TersoffHamannSTM(object):
    """Generate constant height and constant current STM images."""

    def __init__(
        self,
        chg_name="PARCHG",
        min_size=50.0,
        height_tol=2,
        zcut=None,
        use_interpolated=True,
        crop_from_center=False,
        crop_mult=3,
        min_pixel=288,
        interp_step=1.0,
        skew=True,
        extend=0,
        ft_image_path=None,
        ft_image_zoom_factor=5,
        dpi=240,
        cmap="gray",
    ):
        """Initialize class with pathe of PARCHG and other input params."""
        # In original paper, extend used as 1
        chgcar = Chgcar(filename=chg_name)
        self.atoms = chgcar.atoms
        self.dim = chgcar.dim
        self.zmaxp = self.atoms.cart_coords[:, 2].max()
        self.nz = self.dim[2]
        volume = self.atoms.volume
        tmp = chgcar.chg[-1] * volume
        chg = tmp.reshape(self.dim[::-1]).T
        self.chg = chg
        self.cmap = cmap
        self.ft_image_path = ft_image_path
        self.ft_image_zoom_factor = ft_image_zoom_factor
        self.height_tol = height_tol
        self.a = self.atoms.lattice.a
        self.b = self.atoms.lattice.b
        self.c = self.atoms.lattice.c
        self.skew = skew
        self.zcut = zcut
        self.extend = extend
        z_frac_coords = self.atoms.frac_coords[:, 2]
        z_frac_coords_moved = []
        for i in z_frac_coords:
            if i > 0.5:
                i = i - 1
            elif i < -0.5:
                i = i + 1
            z_frac_coords_moved.append(i)
        self.crop_mult = crop_mult
        self.crop_from_center = crop_from_center
        if self.crop_from_center:
            min_size = min_size * self.crop_mult
        self.min_size = min_size
        self.zmaxp = max(np.array(z_frac_coords_moved) * self.c)
        rep_x = int(min_size / self.a) + self.extend
        rep_y = int(min_size / self.b) + self.extend
        self.repeat = [rep_x, rep_y]
        self.scell = self.atoms.make_supercell_matrix([rep_x, rep_y, 1])

        self.use_interpolated = use_interpolated
        self.interp_step = interp_step

        self.min_pixel = min_pixel
        self.dpi = dpi

    def constant_height(
        self,
        filename="testh.png",
    ):
        """Get iso-height image."""
        tol = self.height_tol
        if not self.zcut:
            self.zcut = int((self.zmaxp + tol) / self.c * self.nz)
        # print("zcut", self.zcut, self.repeat)
        info = {}

        img_ext = np.tile(self.chg[:, :, self.zcut], self.repeat)
        if not self.use_interpolated:
            exts = (
                0,
                self.a * self.repeat[0],
                0,
                self.b * (self.repeat[1] - 1),
            )
            plt.close()
            fig, ax = plt.subplots()
            plt.xticks([])
            plt.yticks([])
            if self.skew:
                tmp = 90 - self.atoms.lattice.angles[2]
            else:
                tmp = 0
            data = self.get_plot(
                ax,
                img_ext,
                exts,
                mtransforms.Affine2D().skew_deg(tmp, 0),
            )
            info["data"] = data
            info["img_ext"] = img_ext
            fig.subplots_adjust(bottom=0, top=1, left=0.0, right=1)
            plt.savefig(
                filename, dpi=self.dpi
            )  # , bbox_inches="tight", pad_inches=0.0, dpi=240)
            plt.close()
        else:
            img_ext = self.get_interpolated_data(img_data=img_ext)
            # print ('img_ext',img_ext)
            plt.close()
            # fig, ax = plt.subplots()
            plt.imshow(img_ext, interpolation="none", cmap=self.cmap)
            # ax.set_aspect('equal')
            plt.axis("off")
            plt.savefig(filename, dpi=self.dpi)

            plt.close()

            plt.gca().set_axis_off()
            plt.subplots_adjust(
                top=1, bottom=0, right=1, left=0, hspace=0, wspace=0
            )
            plt.margins(0, 0)
            plt.gca().xaxis.set_major_locator(plt.NullLocator())
            plt.gca().yaxis.set_major_locator(plt.NullLocator())
            # print ('img_ext',img_ext)
            plt.close()

        info["scell"] = self.scell
        info["zcut"] = self.zcut
        info["img_ext"] = img_ext
        if self.crop_from_center:
            im = Image.from_file(filename)

            p_x = int(self.min_pixel / self.crop_mult)
            # print('shape',im.values.shape,p_x)
            im = Image.crop_from_center(
                image_path=filename, target_left=p_x, target_right=p_x
            )
            info["img_ext"] = im
            plt.imshow(im, cmap=self.cmap)
            plt.tight_layout()
            plt.axis("off")
            plt.savefig(filename, dpi=self.dpi)
            plt.close()
        if self.ft_image_path is not None:
            im = Image.from_file(filename).fourier_transform2D(
                zoom_factor=self.ft_image_zoom_factor
            )
            im.save(self.ft_image_path)
        return info

    def constant_current(
        self,
        pc=None,
        ext=0.15,
        filename="testc.png",
        use_interpolated=True,
    ):
        """Return the constant-current cut the charge density."""
        zmax_ind = int(self.zmaxp / self.c * self.nz) + 1
        info = {}
        tol = self.height_tol
        # Find what z value is near the current, and take avergae
        if not self.zcut:
            self.zcut = int((self.zmaxp + tol) / self.c * self.nz)
        zext = int(self.nz * ext)
        zcut_min = self.zcut - zext
        if zcut_min < zmax_ind:
            zcut_min = zmax_ind
        zcut_max = self.zcut + zext
        if pc is None:
            c = np.average(self.chg[:, :, self.zcut])
        else:
            tmp = np.zeros(zcut_max - zcut_min)
            for ii in range(tmp.size):
                tmp[ii] = np.average(self.chg[:, :, zcut_min + ii])
            c = np.linspace(tmp.min(), tmp.max(), 100)[pc]

        # height of iso-current
        img = np.argmin(np.abs(self.chg[:, :, zcut_min:zcut_max] - c), axis=2)
        img_ext = np.tile(img, self.repeat[::-1]) + self.zcut - zext
        if not use_interpolated:
            fig, ax = plt.subplots()
            exts = (
                0,
                self.a * self.repeat[0],
                0,
                self.b * (self.repeat[1] - 1),
            )
            plt.xticks([])
            plt.yticks([])
            if self.skew:
                tmp = 90 - self.atoms.lattice.angles[2]
            else:
                tmp = 0

            data = self.get_plot(
                ax,
                img_ext,
                exts,
                mtransforms.Affine2D().skew_deg(tmp, 0),
            )
            info["data"] = data
            info["img_ext"] = img_ext
            plt.savefig(filename, dpi=self.dpi)
            plt.close()
        else:
            img_ext = self.get_interpolated_data(img_data=img_ext)
            plt.close()
            # fig, ax = plt.subplots()
            plt.imshow(img_ext, interpolation="none", cmap=self.cmap)
            # ax.set_aspect('equal')
            plt.axis("off")
            plt.savefig(filename, dpi=self.dpi)
            plt.close()
        if self.crop_from_center:
            im = Image.from_file(filename)

            p_x = int(self.min_pixel / self.crop_mult)
            # print('shape',im.values.shape,p_x)
            im = Image.crop_from_center(
                image_path=filename, target_left=p_x, target_right=p_x
            )
            plt.imshow(im, cmap=self.cmap)
            plt.axis("off")
            plt.savefig(filename, dpi=self.dpi)
            plt.close()

        info["scell"] = self.scell
        info["zcut"] = self.zcut
        return info

    def get_plot(self, ax, Z, extent, transform):
        """Apply affine transformation."""
        ax.axis("off")
        im = ax.imshow(
            Z,
            interpolation="none",
            origin="lower",
            extent=extent,
            clip_on=True,
            aspect="equal",
        )  # ,cmap=plt.get_cmap('gray')

        trans_data = transform + ax.transData
        im.set_transform(trans_data)

        # display intended extent of the image
        x1, x2, y1, y2 = im.get_extent()
        data = ax.plot(
            [x1, x2, x2, x1, x1],
            [y1, y1, y2, y2, y1],
            "y--",
            transform=trans_data,
        )
        print("min Z and maxZ", np.min(Z), np.max(Z))
        return data

    def get_interpolated_data(self, img_data=[]):
        """Get interpolated data after moving pixels."""
        x = []
        y = []
        zz = []
        xy = []
        step = self.interp_step
        atoms = self.atoms
        for i in range(img_data.shape[0]):
            for j in range(img_data.shape[1]):
                z = img_data[i][j]
                xyz = i * atoms.lattice_mat[0] + j * atoms.lattice_mat[1]
                x.append(xyz[0])
                y.append(xyz[1])
                zz.append(z)
                xy.append([xyz[0], xyz[1]])

        grid_x, grid_y = np.mgrid[
            min(x) : max(x) : step, min(y) : max(y) : step
        ]
        # stepx=(max(x)-min(x))/bins
        # stepy=(max(y)-min(y))/bins
        # grid_x, grid_y = np.mgrid[min(x):max(x):stepx, min(y):max(y):stepy]
        interp = scipy.interpolate.griddata(xy, zz, (grid_x, grid_y)).T
        # plt.scatter(x, y, 1, zz)
        # plt.savefig("ok.png")
        # plt.close()
        return interp


"""
if __name__ == "__main__":
    plt.switch_backend("agg")
    f = "PARCHG"
    t = TersoffHamannSTM(chg_name=f).constant_height()
    t = TersoffHamannSTM(chg_name=f).constant_current()
"""
