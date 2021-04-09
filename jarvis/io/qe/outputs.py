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

    def __init__(self, filename="", data={}, set_key=None):
        """Initialize class."""
        self.filename = filename
        self.data = data
        self.set_key = set_key
        if self.data == {}:
            self.xml_to_dict()

    def xml_to_dict(self):
        """Read XML file."""
        if ".gz" in self.filename:
            import gzip

            with open(self.filename, "rb") as f:
                file_content = f.read()
                self.data = xmltodict.parse(file_content)
        else:
            with open(self.filename) as fd:
                data = xmltodict.parse(fd.read())
                self.data = data
        if self.set_key is None:
            if "step" in self.data["qes:espresso"]:
                self.set_key = "step"
            elif "output" in self.data["qes:espresso"]:
                self.set_key = "output"
            else:
                raise ValueError("Inconsisten QE version.")

    @property
    def final_energy(self):
        """Get final energy."""
        line = self.data["qes:espresso"][self.set_key]
        if isinstance(line, list):
            line = line[-1]
        return float(line["total_energy"]["etot"])

    @property
    def forces(self):
        """Get final forces."""
        line = self.data["qes:espresso"][self.set_key]
        if isinstance(line, list):
            line = line[-1]
        return [
            [float(j) for j in i.split()]
            for i in line["forces"]["#text"].split("\n")
        ]

    @property
    def initial_structure(self):
        """Get input atoms."""
        line = self.data["qes:espresso"]["input"]
        if isinstance(line, list):
            line = line[-1]
        elements = []
        pos = []
        lat = []
        lat.append(
            [float(i) for i in line["atomic_structure"]["cell"]["a1"].split()]
        )
        lat.append(
            [float(i) for i in line["atomic_structure"]["cell"]["a2"].split()]
        )
        lat.append(
            [float(i) for i in line["atomic_structure"]["cell"]["a3"].split()]
        )
        if isinstance(
            line["atomic_structure"]["atomic_positions"]["atom"], list
        ):
            for i in line["atomic_structure"]["atomic_positions"]["atom"]:
                elements.append(i["@name"])
                pos.append([float(j) for j in i["#text"].split()])
            atoms = Atoms(
                elements=elements, coords=pos, lattice_mat=lat, cartesian=True
            )
        else:
            elements = [
                line["atomic_structure"]["atomic_positions"]["atom"]["@name"]
            ]
            pos = [
                [
                    float(j)
                    for j in line["atomic_structure"]["atomic_positions"][
                        "atom"
                    ]["#text"].split()
                ]
            ]
            atoms = Atoms(
                elements=elements, coords=pos, lattice_mat=lat, cartesian=True
            )
        return atoms

    @property
    def final_structure(self):
        """Get final atoms."""
        line = self.data["qes:espresso"][self.set_key]
        if isinstance(line, list):
            line = line[-1]
        elements = []
        pos = []
        lat = []
        lat.append(
            [float(i) for i in line["atomic_structure"]["cell"]["a1"].split()]
        )
        lat.append(
            [float(i) for i in line["atomic_structure"]["cell"]["a2"].split()]
        )
        lat.append(
            [float(i) for i in line["atomic_structure"]["cell"]["a3"].split()]
        )
        if isinstance(
            line["atomic_structure"]["atomic_positions"]["atom"], list
        ):
            for i in line["atomic_structure"]["atomic_positions"]["atom"]:
                elements.append(i["@name"])
                pos.append([float(j) for j in i["#text"].split()])
            atoms = Atoms(
                elements=elements, coords=pos, lattice_mat=lat, cartesian=True
            )
        else:
            elements = [
                line["atomic_structure"]["atomic_positions"]["atom"]["@name"]
            ]
            pos = [
                [
                    float(j)
                    for j in line["atomic_structure"]["atomic_positions"][
                        "atom"
                    ]["#text"].split()
                ]
            ]
            atoms = Atoms(
                elements=elements, coords=pos, lattice_mat=lat, cartesian=True
            )
        return atoms

    @property
    def efermi(self):
        """Get Fermi energy."""
        return float(
            self.data["qes:espresso"]["output"]["band_structure"][
                "fermi_energy"
            ]
        )

    @property
    def functional(self):
        """Get name of DFT functional."""
        return self.data["qes:espresso"]["input"]["dft"]["functional"]

    @property
    def nelec(self):
        """Get number of electrons."""
        return int(
            float(
                self.data["qes:espresso"]["output"]["band_structure"]["nelec"]
            )
        )

    @property
    def indir_gap(self):
        eigs = self.bandstruct_eigvals()
        nelec = self.nelec
        # TODO
        # Check error
        nelec = 3
        return min(eigs[:, nelec + 1]) - max(eigs[:, nelec])

    def bandstruct_eigvals(self, plot=False, filename="band.png"):
        """Get eigenvalues to plot bandstructure."""
        # nbnd = int(
        #    self.data["qes:espresso"]["output"]["band_structure"]["nbnd"]
        # )
        nkpts = int(
            self.data["qes:espresso"]["output"]["band_structure"]["nks"]
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
        # Eigenvalues for each k-point
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
