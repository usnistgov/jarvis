"""Module for analyzing QE outputs."""

from jarvis.core.atoms import Atoms
from collections import OrderedDict
import xmltodict
import numpy as np
import gzip

bohr_to_ang = 0.529177249
hartree_to_ev = 27.2113839


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
            f = gzip.open(self.filename, "rb")
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
        return float(line["total_energy"]["etot"]) * hartree_to_ev

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
    def qe_version(self):
        """Get QE version number."""
        return self.data["qes:espresso"]["general_info"]["creator"]["@VERSION"]

    @property
    def is_spin_orbit(self):
        """Check if spin-orbit coupling in True."""
        tag = self.data["qes:espresso"]["input"]["spin"]["spinorbit"]
        if tag == "true":
            return True
        else:
            return False

    @property
    def is_spin_polarized(self):
        """Check if spin-polarization coupling in True."""
        tag = self.data["qes:espresso"]["output"]["band_structure"]["lsda"]
        if tag == "true":
            return True
        else:
            return False

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
                elements=elements,
                coords=np.array(pos) * bohr_to_ang,
                lattice_mat=np.array(lat) * bohr_to_ang,
                cartesian=True,
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
                elements=elements,
                coords=np.array(pos) * bohr_to_ang,
                lattice_mat=np.array(lat) * bohr_to_ang,
                cartesian=True,
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
                elements=elements,
                coords=np.array(pos) * bohr_to_ang,
                lattice_mat=np.array(lat) * bohr_to_ang,
                cartesian=True,
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
                elements=elements,
                coords=np.array(pos) * bohr_to_ang,
                lattice_mat=np.array(lat) * bohr_to_ang,
                cartesian=True,
            )
        return atoms

    @property
    def efermi(self):
        """Get Fermi energy."""
        return (
            float(
                self.data["qes:espresso"]["output"]["band_structure"][
                    "fermi_energy"
                ]
            )
            * hartree_to_ev
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
    def nkpts(self):
        """Get number of electrons."""
        return int(
            float(self.data["qes:espresso"]["output"]["band_structure"]["nks"])
        )

    @property
    def nbands(self):
        """Get number of bands."""
        return int(
            float(
                self.data["qes:espresso"]["output"]["band_structure"]["nbnd"]
            )
        )

    @property
    def indir_gap(self):
        """Get indirect bandgap."""
        eigs = self.bandstruct_eigvals()  # .T
        nelec = self.nelec
        if not self.is_spin_polarized and nelec % 2 != 0:
            raise ValueError(
                "Odd #electrons cant have band gaps in non-spin-polarized."
            )
        if not self.is_spin_polarized:
            nelec = int(nelec / 2)
        gap = min(eigs[:, nelec]) - max(eigs[:, nelec - 1])
        if gap < 0:
            gap = 0
        return gap

    def bandstruct_eigvals(self, plot=False, filename="band.png"):
        """Get eigenvalues to plot bandstructure."""
        # nbnd = int(
        #    self.data["qes:espresso"]["output"]["band_structure"]["nbnd"]
        # )
        nkpts = self.nkpts
        # int(
        #    self.data["qes:espresso"]["output"]["band_structure"]["nks"]
        # )
        # nbnd = self.nbands
        eigvals = []
        for i in range(nkpts):
            eig = np.array(
                self.data["qes:espresso"]["output"]["band_structure"][
                    "ks_energies"
                ][i]["eigenvalues"]["#text"].split(),
                dtype="float",
            )
            eigvals.append(eig)
        # Eigenvalues for each k-point
        eigvals = np.array(eigvals) * hartree_to_ev
        if plot:
            import matplotlib.pyplot as plt

            for i in eigvals.T:
                plt.plot(i)
            plt.savefig(filename)
            plt.close()
        return eigvals


class ProjHamXml(object):
    """Module to parse projham_K.xml."""

    # Adapted from https://github.com/kfgarrity/TightlyBound.jl
    def __init__(
        self,
        filename="projham_K.xml",
        data=None,
    ):
        """Initialize class."""
        self.filename = filename
        self.data = data
        if data is None:
            self.read()

    def read(self):
        """Read file."""
        if ".gz" in self.filename:
            f = gzip.open(self.filename, "rb")
            file_content = f.read()
            data = xmltodict.parse(file_content)
            self.data = data
        else:
            with open(self.filename) as fd:
                data = xmltodict.parse(fd.read())
                self.data = data

    def get_tight_binding(self):
        """Get tight_binding parameters."""
        tmp_tb = self.data["root"]["tightbinding"]
        nwan = int(float(tmp_tb["nwan"]))
        h1 = np.array(tmp_tb["h1"].split(), dtype="float").reshape(nwan, nwan)
        nonorth = tmp_tb["nonorth"]
        grid = [0, 0, 0]
        if "grid" in list(tmp_tb.keys()):
            grid = np.array(tmp_tb["grid"].split(), dtype="int")
        kweights = np.array(tmp_tb["kweights"].split(), dtype="float")
        nk = int(float(tmp_tb["nk"]))
        kind_arr = np.array(tmp_tb["kind_arr"].split(), dtype="float").reshape(
            3, nk
        )
        hk_lines = tmp_tb["Hk"].split("\n")
        H = np.zeros((nwan, nwan, nk), dtype=complex)
        S = np.zeros((nwan, nwan, nk), dtype=complex)

        for i in hk_lines:
            tmp = i.split()
            if len(tmp) == 7:
                r = int(tmp[0]) - 1
                m = int(tmp[1]) - 1
                n = int(tmp[2]) - 1
                H[m, n, r] = float(tmp[3]) + 1j * float(tmp[4])
                S[m, n, r] = float(tmp[5]) + 1j * float(tmp[6])
        return H, S, h1, kind_arr, kweights, nonorth, grid


# ProjHamXml().get_tight_binding()
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
