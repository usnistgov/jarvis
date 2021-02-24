"""Module to simulate STM with Tershoff-Hamann approach."""
# Reference: https://www.nature.com/articles/s41597-021-00824-y
import matplotlib.pyplot as plt
from jarvis.io.vasp.outputs import Chgcar
import numpy as np
import matplotlib.transforms as mtransforms


class TersoffHamannSTM(object):
    """Generate constant height and constant current STM images."""

    def __init__(self, chg_name="PARCHG", min_size=50.0, skew=True, zcut=None):
        """Initialize class with pathe of PARCHG and other input params."""
        chgcar = Chgcar(filename=chg_name)
        self.atoms = chgcar.atoms
        self.dim = chgcar.dim
        self.zmaxp = self.atoms.cart_coords[:, 2].max()
        self.nz = self.dim[2]
        volume = self.atoms.volume
        tmp = chgcar.chg[-1] * volume
        chg = tmp.reshape(self.dim[::-1]).T
        self.chg = chg
        self.a = self.atoms.lattice.a
        self.b = self.atoms.lattice.b
        self.c = self.atoms.lattice.c
        self.skew = skew
        self.zcut = zcut
        z_frac_coords = self.atoms.frac_coords[:, 2]
        z_frac_coords_moved = []
        for i in z_frac_coords:
            if i > 0.5:
                i = i - 1
            elif i < -0.5:
                i = i + 1
            z_frac_coords_moved.append(i)
        self.zmaxp = max(np.array(z_frac_coords_moved) * self.c)
        rep_x = int(min_size / self.a) + 1
        rep_y = int(min_size / self.b) + 1
        self.repeat = [rep_x, rep_y]
        self.scell = self.atoms.make_supercell_matrix([rep_x, rep_y, 1])

    def constant_height(self, tol=2, filename="testh.png"):
        """Get iso-height image."""
        if not self.zcut:
            self.zcut = int((self.zmaxp + tol) / self.c * self.nz)
        img_ext = np.tile(self.chg[:, :, self.zcut], self.repeat)
        exts = (0, self.a * self.repeat[0], 0, self.b * self.repeat[1])
        plt.close()
        fig, ax = plt.subplots()
        plt.xticks([])
        plt.yticks([])
        if self.skew:
            tmp = 90 - self.atoms.lattice.angles[2]
        else:
            tmp = 0
        data = self.get_plot(
            ax, img_ext, exts, mtransforms.Affine2D().skew_deg(tmp, 0),
        )
        info = {}
        info["img_ext"] = img_ext
        info["data"] = data
        info["scell"] = self.scell
        info["zcut"] = self.zcut
        fig.subplots_adjust(bottom=0, top=1, left=0.0, right=1)
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.0, dpi=240)
        plt.close()

        return info

    def constant_current(self, tol=2, pc=None, ext=0.15, filename="testc.png"):
        """Return the constant-current cut the charge density."""
        zmax_ind = int(self.zmaxp / self.c * self.nz) + 1
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
        fig, ax = plt.subplots()
        exts = (0, self.a * self.repeat[0], 0, self.b * self.repeat[1])
        plt.xticks([])
        plt.yticks([])
        if self.skew:
            tmp = 90 - self.atoms.lattice.angles[2]
        else:
            tmp = 0

        data = self.get_plot(
            ax, img_ext, exts, mtransforms.Affine2D().skew_deg(tmp, 0),
        )
        plt.savefig(filename)
        plt.close()
        info = {}
        info["img_ext"] = img_ext
        info["data"] = data
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


"""
if __name__ == "__main__":
    plt.switch_backend("agg")
    f = "PARCHG"
    t = TersoffHamannSTM(chg_name=f).constant_height()
    t = TersoffHamannSTM(chg_name=f).constant_current()
"""
