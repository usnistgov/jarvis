"""Module for analyzing QE outputs."""


from jarvis.core.atoms import Atoms
from collections import OrderedDict
import xmltodict
import numpy as np


class QEout(object):
    """Module for parsing screen QE output files."""

    def __init__(self, lines=None, filename="qe.out"):
        """Intialize class with filename."""
        self.filename = filename
        if lines is None:
            f = open(filename, "r")
            lns = f.read().splitlines()
            f.close()
            self.lines = lns
        else:
            self.lines = lines

    @classmethod
    def from_dict(self, d={}):
        """Construct from a dictionary."""
        return QEout(lines=d["lines"], filename=d["filename"])

    def to_dict(self):
        """Convert class to a dictionary."""
        d = OrderedDict()
        d["lines"] = self.lines
        d["filename"] = self.filename
        return d

    def get_total_energy(self):
        """Get total energy in Ry."""
        energies = []
        for i in self.lines:
            if "total energy              =" in i:
                print(i)
                energy = float(i.split()[-2])
                energies.append(energy)
        return energies[-1]

    def get_efermi(self):
        """Get fermi energy in eV."""
        efs = []
        for i in self.lines:
            if "the Fermi energy is" in i:
                efs.append(float(i.split()[-2]))
        return efs[-1]

    def get_band_enegies(self):
        """Get band energies in eV."""
        band_energies = []
        for ii, i in enumerate(self.lines):
            if "bands (ev)" in i:
                band_energies.append(
                    [float(j) for j in self.lines[ii + 2].split()]
                )
        return band_energies


class DataFileSchema(object):
    """Module to parse data-file-schema.xml file."""

    def __init__(self, filename="", data={}):
        """Initialize class."""
        self.filename = filename
        self.data = data
        if self.data == {}:
            self.xml_to_dict()

    def xml_to_dict(self):
        """Read XML file."""
        with open(self.filename) as fd:
            data = xmltodict.parse(fd.read())
            self.data = data

    @property
    def final_energy(self):
        """Get final energy."""
        return float(
            self.data["qes:espresso"]["step"][-1]["total_energy"]["etot"]
        )

    @property
    def forces(self):
        """Get final forces."""
        return [
            [float(j) for j in i.split()]
            for i in self.data["qes:espresso"]["step"][-1]["forces"][
                "#text"
            ].split("\n")
        ]

    @property
    def final_structure(self):
        """Get final atoms."""
        elements = []
        pos = []
        lat = []
        lat.append(
            [
                float(i)
                for i in self.data["qes:espresso"]["step"][-1][
                    "atomic_structure"
                ]["cell"]["a1"].split()
            ]
        )
        lat.append(
            [
                float(i)
                for i in self.data["qes:espresso"]["step"][-1][
                    "atomic_structure"
                ]["cell"]["a2"].split()
            ]
        )
        lat.append(
            [
                float(i)
                for i in self.data["qes:espresso"]["step"][-1][
                    "atomic_structure"
                ]["cell"]["a3"].split()
            ]
        )
        for i in self.data["qes:espresso"]["step"][-1]["atomic_structure"][
            "atomic_positions"
        ]["atom"]:
            elements.append(i["@name"])
            pos.append([float(j) for j in i["#text"].split()])
        atoms = Atoms(
            elements=elements, coords=pos, lattice_mat=lat, cartesian=True
        )
        return atoms

    def bandstruct_eigvals(self, plot=False, filename="band.png"):
        """Get eigenvalues to plot bandstructure."""
        # nbnd = int(
        #    self.data["qes:espresso"]["output"]["band_structure"]["nbnd"]
        # )
        nkpts = int(
            self.data["qes:espresso"]["output"]["band_structure"][
                "starting_k_points"
            ]["nk"]
        )
        eigvals = []
        for i in range(nkpts):
            eig = [
                float(j)
                for j in self.data["qes:espresso"]["output"]["band_structure"][
                    "ks_energies"
                ][i]["eigenvalues"]["#text"]
                .split("\n")[0]
                .split()
            ]
            eigvals.append(eig)
        eigvals = np.array(eigvals)
        if plot:
            import matplotlib.pyplot as plt

            for i in eigvals.T:
                plt.plot(i)
            plt.savefig(filename)
            plt.close()
        return eigvals


"""
if __name__ == "__main__":
    en = QEout("qe.out").get_total_energy()
    print((en))
    assert en == -19.11812163

    en = QEout("qe.out").get_band_enegies()
    print((en), len(en))
    assert en[0][0] == -5.8325
    en = QEout("qe.out").get_efermi()
    print((en))
    assert en == 6.4236
"""
