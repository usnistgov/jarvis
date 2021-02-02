"""Class for analyzing  WT.out file."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.switch_backend("agg")


class WTOut(object):
    """Construct WT.out related object."""

    def __init__(self, path=""):
        """Initialize with path to WT.out file."""
        self.path = path

    def get_z2_index(self):
        """Get Z2 index."""
        f = open(self.path, "r")
        lines = f.read().splitlines()
        f.close()
        strng = 0
        for ii, i in enumerate(lines):
            if "# z2 number for 6 planes" in i:
                for j in range(6):
                    tmp = lines[ii + j + 1].split()

                    if (
                        tmp[0] == "k1=0.0,"
                        or tmp[0] == "k2=0.0,"
                        or tmp[0] == "k3=0.0,"
                        or tmp[0] == "k3=0.5,"
                        or tmp[0] == "k2=0.5,"
                        or tmp[0] == "k1=0.5,"
                    ):
                        val = tmp[3]
                        strng = strng + float(val)
                    if tmp[0] == "k1=0.5,":
                        weak1 = tmp[3]
                    if tmp[0] == "k2=0.5,":
                        weak2 = tmp[3]
                    if tmp[0] == "k3=0.5,":
                        weak3 = tmp[3]
        index = (
            str(int(strng % 2))
            + str(";")
            + str(weak1)
            + str(",")
            + str(weak2)
            + str(",")
            + str(weak3)
        )
        return index

    def get_chern_number(self):
        """Get Chern index."""
        f = open(self.path, "r")
        lines = f.read().splitlines()
        f.close()
        chrn = []
        try:
            for j, jj in enumerate(lines):
                if "# Chern number for 6 planes" in jj:
                    for k in range(1, 7):
                        tmp = float(lines[j + k].split(":")[1])
                        if tmp not in chrn:
                            chrn.append(tmp)
        except Exception:
            pass
        return chrn


def parse_chern_dat(
    chern_dat="wanniercenter3D_Chern.dat", filename="mychern.png"
):
    """Plot wanniercenter3D_Chern.dat file."""
    x = np.loadtxt(chern_dat)
    if filename is not None:
        the_grid = GridSpec(3, 2)
        plt.rcParams.update({"font.size": 18})
        plt.figure(figsize=(12, 12))

        plt.subplot(the_grid[0, 0])
        plt.title("(a) k$_1$=0.0")
        plt.xlabel("k$_2$")
        plt.plot(x[:, 0], x[:, 1], ".")
        plt.subplot(the_grid[0, 1])
        plt.title("(b) k$_1$=0.5")
        plt.xlabel("k$_2$")
        plt.plot(x[:, 0], x[:, 2], ".")

        plt.subplot(the_grid[1, 0])
        plt.xlabel("k$_1$")
        plt.title("(c) k$_2$=0.0")
        plt.plot(x[:, 0], x[:, 3], ".")
        plt.subplot(the_grid[1, 1])
        plt.xlabel("k$_1$")
        plt.title("(d) k$_2$=0.5")
        plt.plot(x[:, 0], x[:, 4], ".")

        plt.subplot(the_grid[2, 0])
        plt.title("(e) k$_3$=0.0")
        plt.xlabel("k$_2$")
        plt.plot(x[:, 0], x[:, 5], ".")
        plt.subplot(the_grid[2, 1])
        plt.title("(f) k$_3$=0.5")
        plt.xlabel("k$_2$")
        plt.plot(x[:, 0], x[:, 6], ".")
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()
    return x


def parse_nodes_dat(fname="Nodes.dat"):
    """Parse Nodedat file."""
    f = open(fname, "r")
    lines = f.read().splitlines()
    f.close()
    nodes = []
    for i in lines:
        if "#" not in i:
            nodes.append(i.split())
    nodes = np.array(nodes, dtype="float")
    return nodes


"""
if __name__ == "__main__":
    wt = "WT.out"
    z2 = WTOut(path=wt).get_z2_index()
    print(z2)
    chrn = WTOut(path=wt).get_chern_number()
    print(chrn)
"""
