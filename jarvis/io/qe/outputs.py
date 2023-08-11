"""Module for analyzing QE outputs."""

from jarvis.core.atoms import Atoms
from collections import OrderedDict
import xmltodict
import numpy as np
import gzip
import scipy.linalg as la
import copy

bohr_to_ang = 0.529177249
hartree_to_ev = 27.2113839
ryd_to_ev = hartree_to_ev / 2.0


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
                energy = float(i.split()[-2])
                energies.append(energy)
        return float(energies[-1]) * ryd_to_ev

    def get_efermi(self):
        """Get fermi energy in eV."""
        efs = []
        for i in self.lines:
            if "the Fermi energy is" in i:
                efs.append(float(i.split()[-2]))
        return efs[-1]

    @property
    def job_done(self):
        """Check if job is completed."""
        done = False
        for i in self.lines:
            if "JOB DONE." in i:
                done = True
        return done

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
            if "output" in self.data["qes:espresso"]:
                self.set_key = "output"
            elif "step" in self.data["qes:espresso"]:
                self.set_key = "step"
            else:
                raise ValueError("Inconsisten QE version.")

    @property
    def final_energy(self):
        """Get final energy."""
        # Energy in eV already ??
        line = self.data["qes:espresso"][self.set_key]
        if isinstance(line, list):
            line = line[-1]
        return float(line["total_energy"]["etot"])  # * hartree_to_ev

    @property
    def num_atoms(self):
        """Get total number of atoms."""
        return self.final_structure.num_atoms

    @property
    def energy_per_atom(self):
        """Get final energy per atom."""
        return self.final_energy / self.num_atoms

    @property
    def final_energy_breakdown(self):
        """Get final energy."""
        line = self.data["qes:espresso"][self.set_key]
        if isinstance(line, list):
            line = line[-1]
        tmp = line["total_energy"]
        for i, j in tmp.items():
            tmp[i] = float(j) * hartree_to_ev
        return tmp

    @property
    def forces(self):
        """Get final forces."""
        line = self.data["qes:espresso"][self.set_key]
        if isinstance(line, list):
            line = line[-1]
        return np.array(
            [
                [float(j) for j in i.split()]
                for i in line["forces"]["#text"].split("\n")
            ]
        ) * (hartree_to_ev / bohr_to_ang)

    @property
    def stress(self):
        """Get final stress."""
        line = self.data["qes:espresso"][self.set_key]
        if isinstance(line, list):
            line = line[-1]
        return np.array(
            [
                [float(j) for j in i.split()]
                for i in line["stress"]["#text"].split("\n")
            ]
        ) * (hartree_to_ev / bohr_to_ang**3)

    @property
    def magnetization(self):
        """Get final magnetization."""
        line = self.data["qes:espresso"][self.set_key]
        if isinstance(line, list):
            line = line[-1]
        return float(line["magnetization"]["total"])

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
        if self.is_spin_polarized:
            return [
                int(
                    float(
                        self.data["qes:espresso"]["output"]["band_structure"][
                            "nbnd_up"
                        ]
                    )
                ),
                int(
                    float(
                        self.data["qes:espresso"]["output"]["band_structure"][
                            "nbnd_dw"
                        ]
                    )
                ),
            ]
        else:
            return int(
                float(
                    self.data["qes:espresso"]["output"]["band_structure"][
                        "nbnd"
                    ]
                )
            )

    @property
    def indir_gap(self):
        """Get indirect bandgap."""
        eigs = self.bandstruct_eigvals()  # .T
        nelec = self.nelec
        if not self.is_spin_polarized:
            nelec = int(nelec / 2)
        gap = min(eigs[:, nelec]) - max(eigs[:, nelec - 1])
        if not self.is_spin_polarized and nelec % 2 != 0 and gap > 0.1:
            raise ValueError(
                "Odd #electrons cant have band gaps in non-spin-polarized.",
                self.is_spin_polarized,
                nelec,
            )
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

    def dos(self, smearing=0.2, plot=False, filename="dos.png"):
        """Density of states."""
        """Based on sum of gaussians with smearing as given"""

        # TODO: make work nicely for spin-polarized case,
        # with minority spins plotted negative.

        nkpts = self.nkpts
        eigvals = []
        kweight = []
        for i in range(nkpts):
            eig = np.array(
                self.data["qes:espresso"]["output"]["band_structure"][
                    "ks_energies"
                ][i]["eigenvalues"]["#text"].split(),
                dtype="float",
            )
            eigvals.append(eig)
            kweight.append(
                float(
                    self.data["qes:espresso"]["output"]["band_structure"][
                        "ks_energies"
                    ][i]["k_point"]["@weight"]
                )
            )

        efermi = self.efermi
        eigvals = np.array(eigvals) * hartree_to_ev - efermi
        kweight = np.array(kweight)

        minval = np.min(np.array(eigvals))
        maxval = np.max(np.array(eigvals))

        energies = np.arange(minval - 0.5, maxval + 0.5, 0.01)
        # de = 0.01
        norm = (1 / 2.0 / np.pi / smearing**2) ** 0.5
        DOS = np.zeros(np.shape(energies)[0])

        for k in range(nkpts):
            for e in eigvals[k, :]:
                DOS += (
                    kweight[k]
                    * norm
                    * np.exp(-0.5 * (energies - e) ** 2 / smearing**2)
                )

        # print("k ", np.sum(kweight))
        # print("DOS ", np.sum(DOS*de),
        # " should be close to ", self.nbands * np.sum(kweight))

        if plot:
            import matplotlib.pyplot as plt

            plt.plot(energies, DOS)
            plotmin = max(-10.0, minval)

            plt.plot(
                [0.0, 0.0],
                [0.0, np.max(DOS) * 1.1],
                color="black",
                linestyle="dashed",
            )
            plt.ylim([0, np.max(DOS) * 1.1])
            plt.xlim([plotmin, maxval])
            plt.ylabel("DOS (eV$^{-1}$)")
            plt.xlabel("Energy - E$_F$ (eV)")
            #    plt.show()
            plt.savefig(filename)
            plt.close()
        return energies, DOS


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

        (
            H,
            S,
            h1,
            kind_arr,
            kweights,
            nonorth,
            grid,
            scf,
            nelec,
        ) = self.get_tight_binding()
        self.H = H
        self.S = S
        self.h1 = h1
        self.kind_arr = kind_arr
        self.kweights = kweights
        self.nonorth = nonorth
        self.grid = grid
        self.nk = np.shape(kind_arr)[0]
        self.nwan = np.shape(H)[1]
        self.scf = scf
        self.nelec = nelec

        A, coords, types, nat = self.get_crystal()
        self.A = A
        self.coords = coords
        self.types = types
        self.nat = nat

        scf = self.data["root"]["scf"]
        if scf is True:
            print("error scf ", scf)

        # info for deciding which orbital goes with which index.
        self.atomdata = dict()
        self.atomdata["H"] = ["s"]
        self.atomdata["Li"] = ["s", "p"]
        self.atomdata["Be"] = ["s", "p"]
        self.atomdata["B"] = ["s", "p"]
        self.atomdata["C"] = ["s", "p"]
        self.atomdata["N"] = ["s", "p"]
        self.atomdata["O"] = ["s", "p"]
        self.atomdata["F"] = ["s", "p"]

        self.atomdata["Na"] = ["s", "p"]
        self.atomdata["Mg"] = ["s", "p"]
        self.atomdata["Al"] = ["s", "p"]

        self.atomdata["Si"] = ["s", "p"]
        self.atomdata["P"] = ["s", "p"]
        self.atomdata["S"] = ["s", "p"]
        self.atomdata["Cl"] = ["s", "p"]
        self.atomdata["K"] = ["s", "d", "p"]
        self.atomdata["Ca"] = ["s", "d", "p"]
        self.atomdata["Sc"] = ["s", "d", "p"]
        self.atomdata["Ti"] = ["s", "d", "p"]
        self.atomdata["V"] = ["s", "d", "p"]
        self.atomdata["Cr"] = ["s", "d", "p"]
        self.atomdata["Mn"] = ["s", "d", "p"]
        self.atomdata["Fe"] = ["s", "d", "p"]
        self.atomdata["Co"] = ["s", "d", "p"]
        self.atomdata["Ni"] = ["s", "d", "p"]
        self.atomdata["Cu"] = ["s", "d", "p"]
        self.atomdata["Zn"] = ["s", "d", "p"]
        self.atomdata["Ga"] = ["s", "d", "p"]
        self.atomdata["Ge"] = ["s", "p"]
        self.atomdata["As"] = ["s", "p"]
        self.atomdata["Se"] = ["s", "p"]
        self.atomdata["Br"] = ["s", "p"]

        self.atomdata["Rb"] = ["s", "d", "p"]
        self.atomdata["Sr"] = ["s", "d", "p"]
        self.atomdata["Y"] = ["s", "d", "p"]
        self.atomdata["Zr"] = ["s", "d", "p"]
        self.atomdata["Nb"] = ["s", "d", "p"]
        self.atomdata["Mo"] = ["s", "d", "p"]
        self.atomdata["Tc"] = ["s", "d", "p"]
        self.atomdata["Ru"] = ["s", "d", "p"]
        self.atomdata["Rh"] = ["s", "d", "p"]
        self.atomdata["Pd"] = ["s", "d", "p"]
        self.atomdata["Ag"] = ["s", "d", "p"]
        self.atomdata["Cd"] = ["s", "d", "p"]
        self.atomdata["In"] = ["s", "d", "p"]
        self.atomdata["Sn"] = ["s", "p"]
        self.atomdata["Sb"] = ["s", "p"]
        self.atomdata["Te"] = ["s", "p"]
        self.atomdata["I"] = ["s", "p"]
        self.atomdata["Cs"] = ["s", "d", "p"]
        self.atomdata["Ba"] = ["s", "d", "p"]
        self.atomdata["La"] = ["s", "d"]
        self.atomdata["Hf"] = ["s", "d", "p"]
        self.atomdata["Ta"] = ["s", "d", "p"]
        self.atomdata["W"] = ["s", "d", "p"]
        self.atomdata["Re"] = ["s", "d", "p"]
        self.atomdata["Os"] = ["s", "d", "p"]
        self.atomdata["Ir"] = ["s", "d", "p"]
        self.atomdata["Pt"] = ["s", "d", "p"]
        self.atomdata["Au"] = ["s", "d", "p"]
        self.atomdata["Hg"] = ["s", "d", "p"]
        self.atomdata["Tl"] = ["s", "d", "p"]
        self.atomdata["Pb"] = ["s", "p"]
        self.atomdata["Bi"] = ["s", "p"]

    def get_crystal(self):
        """Get crystal info."""
        tmp_tb = self.data["root"]["crystal"]
        A = np.reshape(np.array(tmp_tb["A"].split(), dtype="float"), (3, 3))
        nat = int(float(tmp_tb["nat"]))
        coords = np.reshape(
            np.array(tmp_tb["coords"].split(), dtype="float"), (nat, 3)
        )
        types = tmp_tb["types"].split()

        if len(types) != nat:
            print("error loading crystal ", nat, types)
        if np.shape(coords)[0] != nat:
            print("error loading crystal ", nat, np.shape(coords)[0])

        return A, coords, types, nat

    def get_tight_binding(self):
        """Get tight_binding parameters."""
        t = self.data["root"]["scf"]
        if t == "false":
            scf = False
        else:
            scf = True

        nelec = float(self.data["root"]["nelec"])

        tmp_tb = self.data["root"]["tightbinding"]
        nwan = int(float(tmp_tb["nwan"]))

        if "h1" in tmp_tb:
            h1 = np.array(tmp_tb["h1"].split(), dtype="float").reshape(
                nwan, nwan
            )
        else:
            h1 = np.zeros((nwan, nwan), dtype=float)

        t = tmp_tb["nonorth"]
        if t == "false":
            nonorth = False
        else:
            nonorth = True

        grid = [0, 0, 0]
        if "grid" in list(tmp_tb.keys()):
            grid = np.array(tmp_tb["grid"].split(), dtype="int")

        kweights = np.array(tmp_tb["kweights"].split(), dtype="float")
        nk = int(float(tmp_tb["nk"]))
        kind_arr = np.array(tmp_tb["kind_arr"].split(), dtype="float").reshape(
            nk, 3
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
        return H, S, h1, kind_arr, kweights, nonorth, grid, scf, nelec

    def calculate_eigenvalues(self, kpoint=None, kind=-1):
        """Calculate eigenvalues."""
        if self.scf is True:
            print("warning, not accurate for scf=True")

        if kind == -1 and kpoint is not None:  # must find kpoint
            for k in range(self.nk):
                if np.sum(np.abs(kpoint - self.kind_arr[k, :])) < 1e-3:
                    kind = k
                    break
            if kind == -1:
                print("warning, kpoint not found ", kpoint)

        elif kind == -1 and kpoint is None:
            kind = 0
            print("warning, no k-point specified, setting kind to 0")

            # else kind is correctly specified

        if kind >= self.nk:
            print("warning, kind too large ", kind)

        hk = self.H[:, :, kind]
        hk = 0.5 * (hk + np.conj(hk.T))
        if self.nonorth is False:
            vals, vects = la.eigh(hk)
            return vals, vects
        elif self.nonorth is True:
            sk = self.S[:, :, kind]
            sk = 0.5 * (sk + np.conj(sk.T))
            vals, vects = la.eigh(hk, b=sk)
            return vals, vects

    def solve_ham(self, proj=None):
        """Solve hamiltonian."""
        VALS = np.zeros((self.nk, self.nwan), dtype=float)

        if proj is not None:
            nproj = len(proj)
            PROJ = np.zeros((self.nk, self.nwan, nproj), dtype=float)
        else:
            PROJ = None
            nproj = 0

        for k in range(self.nk):
            vals, vects = self.calculate_eigenvalues(kind=k)
            VALS[k, :] = vals

            sk = self.S[:, :, k]
            if proj is not None:
                for ip, p in enumerate(proj):
                    for pind in p:
                        for a in range(self.nwan):
                            for b in range(self.nwan):
                                t = vects[pind, a] * np.conj(vects[b, a])
                                PROJ[k, a, ip] += 0.5 * np.real(
                                    t * sk[b, pind]
                                    + np.conj(t) * np.conj(sk[b, pind])
                                )

        return VALS, PROJ

    # figure our correspondence between orbitals and indicies.
    def count_orbs(self):
        """Count orbitals."""
        ORBS = []
        c = 0
        for a in range(self.nat):
            t = self.types[a]
            orbs = self.atomdata[t]
            for o in orbs:
                if o == "s":
                    ORBS.append([c, t, o])
                    c += 1
                elif o == "p":
                    ORBS.append([c, t, o])
                    ORBS.append([c + 1, t, o])
                    ORBS.append([c + 2, t, o])
                    c += 3
                elif o == "d":
                    ORBS.append([c, t, o])
                    ORBS.append([c + 1, t, o])
                    ORBS.append([c + 2, t, o])
                    ORBS.append([c + 3, t, o])
                    ORBS.append([c + 4, t, o])
                    c += 5
                elif o == "f":
                    ORBS.append([c, t, o])
                    ORBS.append([c + 1, t, o])
                    ORBS.append([c + 2, t, o])
                    ORBS.append([c + 3, t, o])
                    ORBS.append([c + 4, t, o])
                    ORBS.append([c + 5, t, o])
                    ORBS.append([c + 6, t, o])
                    c += 7
        return ORBS

    # figure our orbitials to project onto from inputs
    def decide_projection(self, proj_atoms=None, proj_orbs=None):
        """Decide projections."""
        ORBS = self.count_orbs()

        if proj_atoms is None and proj_orbs is None:
            ntypes = len(set(self.types))
            if ntypes > 1:
                t = list(set(self.types))
                names = t
                proj_atoms = []
                proj_orbs = []
                for tt in t:
                    proj_atoms.append([tt])
                    proj_orbs.append(["s", "p", "d"])
            else:
                names = ["s", "p", "d"]
                proj_orbs = [["s"], ["p"], ["d"]]
                t = list(set(self.types))
                proj_atoms = [t, t, t]
        else:
            if proj_orbs is None:
                names = copy.copy(proj_atoms)
                proj_atoms = []
                proj_orbs = []
                for tt in names:
                    proj_atoms.append([tt])
                    proj_orbs.append(["s", "p", "d"])
            elif proj_atoms is None:
                names = copy.copy(proj_orbs)
                proj_orbs = []
                proj_atoms = []
                t = list(set(self.types))
                for tt in names:
                    proj_orbs.append([tt])
                    proj_atoms.append(t)

            else:
                names = []
                proj_atoms_t = copy.copy(proj_atoms)
                proj_orbs_t = copy.copy(proj_orbs)

                proj_atoms = []
                proj_orbs = []
                for a, o in zip(proj_atoms_t, proj_orbs_t):
                    names.append(a + "_" + o)
                    proj_atoms.append([a])
                    proj_orbs.append([o])

                    #        print("proj_atoms ", proj_atoms)
                    #        print("proj_orbs ", proj_orbs)

        num_proj = len(proj_atoms)
        proj = []
        for cp in range(num_proj):
            proj.append([])

            for c, o in enumerate(ORBS):
                if o[1] in proj_atoms[cp] and o[2] in proj_orbs[cp]:
                    proj[-1].append(c)

        return proj, names

    def dos(
        self,
        smearing=0.3,
        npts=500,
        proj_atoms=None,
        proj_orbs=None,
        do_proj=True,
    ):
        """Get DOS."""
        if do_proj is False:
            proj = None
            names = None
        else:
            proj, names = self.decide_projection(
                proj_atoms=proj_atoms, proj_orbs=proj_orbs
            )

        VALS, PROJ = self.solve_ham(proj=proj)

        VALS = VALS * ryd_to_ev

        vmin = np.min(VALS) - smearing * 5
        vmax = np.max(VALS) + smearing * 5

        de = vmax - vmin

        energies = np.arange(vmin, vmax, de / npts)
        npts = len(energies)

        dos = np.zeros(npts)

        W = np.tile(self.kweights, (self.nwan, 1)).T / 2.0

        for c, e in enumerate(energies):
            dos[c] = np.sum(
                np.exp(-0.5 * (VALS[:, :] - e) ** 2 / smearing**2) * W
            )

        dos = dos / smearing / (2.0 * np.pi) ** 0.5

        if do_proj:
            nproj = len(names)
            pdos = np.zeros((npts, nproj))

            for i in range(nproj):
                for c, e in enumerate(energies):
                    pdos[c, i] = np.sum(
                        PROJ[:, :, i]
                        * np.exp(-0.5 * (VALS[:, :] - e) ** 2 / smearing**2)
                        * W
                    )

            pdos = pdos / smearing / (2.0 * np.pi) ** 0.5

        de = energies[1] - energies[0]
        print("Int DOS = ", np.sum(dos) * de)

        if do_proj:
            for i in range(nproj):
                print(names[i], " Int pDOS = ", np.sum(pdos[:, i]) * de)

        occ = np.cumsum(dos * de)

        for i in range(npts):
            if occ[i] > self.nelec / 2.0:
                fermi_ind = i
                break

        print(
            "Int occupied DOS (only exact as npts goes to inf) = ",
            np.sum(dos[0 : fermi_ind + 1]) * de * 2.0,
        )

        energies = energies - energies[fermi_ind]  # shift fermi energy to zero

        return energies, dos, pdos, names


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

# p = ProjHamXml("/home/kfg/projham_K.xml.gz")
# print("A")
# print(p.A)
# print("coords")
# print(p.coords)
# print("types")
# print(p.types)

# energies, dos, pdos, names = p.dos()

# import matplotlib.pyplot as plt

# plt.plot(energies, dos, "b")
# plt.plot(energies, pdos[:,0], "r")
# plt.plot(energies, pdos[:,1], "g")
# plt.show()
