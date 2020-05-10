"""
Modules for analzing VASP outputs
"""

from scipy.constants import physical_constants, speed_of_light
from jarvis.io.vasp.vasp_constant import *
from jarvis.core.atoms import Atoms
import numpy as np
import xmltodict
from collections import OrderedDict
from jarvis.core.atoms import Atoms
from jarvis.core.specie import Specie
import numpy as np
import xmltodict
from collections import OrderedDict
from jarvis.core.kpoints import Kpoints3D as Kpoints
import sys
from jarvis.io.vasp.inputs import Poscar
import matplotlib
from matplotlib import pyplot as plt


class Chgcar(object):
    def __init__(
        self,
        filename="",
        atoms=None,
        chg=[],
        chgdif=None,
        aug=None,
        augdiff=None,
        nsets=1,
    ):
        """
        Class handling VASP CHGCAR file data
        """
        self.filename = filename
        self.atoms = atoms
        self.chg = chg
        self.chgdif = chgdif
        self.aug = aug
        self.augdiff = augdiff
        self.nsets = nsets
        if self.atoms is None:
            self.read_file()

    def is_spin_polarized(self):
        if self.nsets == 2:
            return True
        else:
            return False

    def is_spin_orbit(self):
        if self.nsets == 4:
            return True
        else:
            return False

    def read_file(self):
        f = open(self.filename, "r")
        lines = f.read()
        f.close()
        self.atoms = Poscar.from_string(lines).atoms
        volume = self.atoms.volume
        text = lines.splitlines()
        for ii, i in enumerate(text):
            if i == "":
                ng_line = text[ii + 1]
                ng = [int(j) for j in text[ii + 1].split()]
                found = ii
        nsets = 0
        for i in text:
            if "augmentation occupancies   1 " in i:
                nsets = nsets + 1
        self.nsets = nsets
        if self.is_spin_orbit():
            ValueError("Not implemeted for spin-orbit calculations yet")
        # print ('nsets=',nsets)

        start = found + 2
        ngs = int(ng[0] * ng[1] * ng[2])
        if ngs % 5 == 0:
            nlines = int(ngs / 5.0)
        else:
            nlines = int(ngs / 5.0) + 1
        end = nlines + start  # +1

        for ii, i in enumerate(text):
            if text[ii] == ng_line:
                start = ii + 1
                end = start + nlines
                chg = self.chg_set(text, start, end, volume, ng)
                self.chg.append(chg)

    def chg_set(self, text, start, end, volume, ng):
        chg = np.empty(ng)
        # print ('start,end',start,end)
        # if 'augmentation occupancies' not in text[end]:
        #   raise ValueError('nlines is set incoorect',text[end])

        lines_0 = text[start:end]
        tmp = []
        for i in lines_0:
            for j in i.split():
                if j != "":
                    tmp.append(float(j))
        # print ('tmp=',len(tmp),tmp[0],tmp[-1])
        # print ('lines_0',len(lines_0),ng,ng[0]*ng[1]*ng[2],lines_0[0].split())
        tmp = np.array(tmp).reshape(ng)
        tmp = tmp / volume
        return tmp


class Vasprun(object):
    def __init__(self, filename="vasprun.xml", data={}):
        """
        Class handling VASP vasprun.xml file data
        """
        self._filename = filename
        self._data = data
        self.ionic_steps = None
        self.electronic_steps = None
        self.input_parameters = None
        if self._data == {}:
            self.xml_to_dict()

    def xml_to_dict(self):
        with open(self._filename) as fd:
            data = xmltodict.parse(fd.read())
            self._data = data
            self.ionic_steps = data["modeling"]["calculation"]
            if type(self.ionic_steps) is not list:
                self.ionic_steps = [self.ionic_steps]
            if self.input_parameters is None:
                self.input_parameters = self.all_input_parameters

    @property
    def final_energy(self):
        return float(self.ionic_steps[-1]["scstep"][-1]["energy"]["i"][11]["#text"])

    @property
    def efermi(self):
        return float(self.ionic_steps[-1]["dos"]["i"]["#text"])

    @property
    def num_atoms(self):
        return self._data["modeling"]["atominfo"]["atoms"]

    @property
    def num_types(self):
        return int(self._data["modeling"]["atominfo"]["types"])

    @property
    def dielectric_loptics(self):
        tmp = self.ionic_steps[-1]["dielectricfunction"]["real"]["array"]["set"]["r"]
        reals = []
        for i in range(len(tmp)):
            reals.append([float(j) for j in tmp[i].split()])

        tmp = self.ionic_steps[-1]["dielectricfunction"]["imag"]["array"]["set"]["r"]
        imags = []
        for i in range(len(tmp)):
            imags.append([float(j) for j in tmp[i].split()])
        reals = np.array(reals)
        imags = np.array(imags)
        return reals, imags

    @property
    def avg_absorption_coefficient(self, max_axis=3):
        eV_to_recip_cm = 1.0 / (
            physical_constants["Planck constant in eV s"][0] * speed_of_light * 1e2
        )
        real, imag = self.dielectric_loptics
        energies = real[:, 0]
        epsilon_1 = np.mean(real[:, 1:max_axis], axis=1)
        epsilon_2 = np.mean(imag[:, 1:max_axis], axis=1)
        absorption = (
            2
            * np.pi
            * np.sqrt(2.0)
            * eV_to_recip_cm
            * energies
            * np.sqrt(-epsilon_1 + np.sqrt(epsilon_1 ** 2 + epsilon_2 ** 2))
        )
        return energies, absorption

    @property
    def dfpt_data(self, fc_mass=True):
        info = {}
        hessian = []
        data = self._data
        for i in (data["modeling"]["calculation"]["dynmat"]["varray"])[0]["v"]:
            hessian.append(i.split())
        hessian = np.array(hessian, dtype="double")
        struct = self.all_structures[-1]
        natoms = struct.num_atoms
        force_constants = np.zeros((natoms, natoms, 3, 3), dtype="double")
        for i in range(natoms):
            for j in range(natoms):
                force_constants[i, j] = hessian[
                    i * 3 : (i + 1) * 3, j * 3 : (j + 1) * 3
                ]
        masses = [Specie(i).atomic_mass for i in struct.elements]
        if fc_mass == True:
            for i in range(natoms):
                for j in range(natoms):
                    force_constants[i, j] *= -np.sqrt(masses[i] * masses[j])
        born_charges = []
        for n in range(natoms):
            born_charges.append(
                [
                    i.split()
                    for i in (data["modeling"]["calculation"]["array"]["set"])[n]["v"]
                ]
            )
        born_charges = np.array(born_charges, dtype="double")
        phonon_eigenvals = np.array(
            data["modeling"]["calculation"]["dynmat"]["v"]["#text"].split(),
            dtype="double",
        )
        eigvecs = np.array(
            [
                i.split()
                for i in (data["modeling"]["calculation"]["dynmat"]["varray"][1]["v"])
            ],
            dtype="float",
        )
        phonon_eigenvectors = []
        for ev in eigvecs:
            phonon_eigenvectors.append(np.array(ev).reshape(natoms, 3))
        phonon_eigenvectors = np.array(phonon_eigenvectors, dtype="float")
        epsilon = {}
        for i in data["modeling"]["calculation"]["varray"]:
            if "epsilon" in i["@name"]:
                epsilon[i["@name"]] = np.array(
                    [j.split() for j in i["v"]], dtype="float"
                )
        info["born_charges"] = born_charges
        info["phonon_eigenvectors"] = phonon_eigenvectors
        info["epsilon"] = epsilon
        info["phonon_eigenvalues"] = phonon_eigenvals
        info["masses"] = masses
        return info

    @property
    def get_dir_gap(self):
        if not self.is_spin_polarized:
            spin_channels = 2
            up_gap = "na"
            dn_gap = "na"
            if self.is_spin_orbit:
                spin_channels = 1
            levels = int(
                float(self.all_input_parameters["NELECT"]) / float(spin_channels)
            )
            ups = self.eigenvalues[0][:, :, 0][:, levels]
            dns = self.eigenvalues[0][:, :, 0][:, levels - 1]
            gap = min(ups - dns)
        else:
            tmp = np.concatenate(
                (self.eigenvalues[0][:, :, 0], self.eigenvalues[1][:, :, 0]), axis=1
            )
            cat = np.sort(tmp, axis=1)
            nelect = int(float(self.all_input_parameters["NELECT"]))
            ups = cat[:, nelect]
            dns = cat[:, nelect - 1]
            gap = min(ups - dns)

        return gap

    @property
    def get_indir_gap(self):
        if not self.is_spin_polarized:
            spin_channels = 2
            up_gap = "na"
            dn_gap = "na"
            if self.is_spin_orbit:
                spin_channels = 1
            levels = int(
                float(self.all_input_parameters["NELECT"]) / float(spin_channels)
            )
            print("levels", levels)
            gap = min(self.eigenvalues[0][:, :, 0][:, levels]) - max(
                self.eigenvalues[0][:, :, 0][:, levels - 1]
            )

        if self.is_spin_polarized:
            print(self.eigenvalues[0][:, :, 0], self.eigenvalues[0][:, :, 0].shape)
            print(self.eigenvalues[1][:, :, 0], self.eigenvalues[1][:, :, 0].shape)
            tmp = np.concatenate(
                (self.eigenvalues[0][:, :, 0], self.eigenvalues[1][:, :, 0]), axis=1
            )
            cat = np.sort(tmp, axis=1)
            nelect = int(float(self.all_input_parameters["NELECT"]))
            gap = min(cat[:, nelect]) - max(cat[:, nelect - 1])
        return gap

    @property
    def elements(self):
        elements = [
            self._data["modeling"]["atominfo"]["array"][0]["set"]["rc"][i]["c"][0]
            for i in range(
                len(self._data["modeling"]["atominfo"]["array"][0]["set"]["rc"])
            )
        ]
        if len(elements) != self.num_atoms:
            ValueError("Number of atoms is  not equal to number of elements")
        elements = [str(i) for i in elements]
        return elements

    def vrun_structure_to_atoms(self, s={}):
        lattice_mat = np.array(
            [[float(j) for j in i.split()] for i in s["crystal"]["varray"][0]["v"]]
        )
        frac_coords = np.array(
            [[float(j) for j in i.split()] for i in s["varray"]["v"]]
        )
        elements = self.elements
        atoms = Atoms(
            lattice_mat=lattice_mat,
            elements=elements,
            coords=frac_coords,
            cartesian=False,
        )
        return atoms

    @property
    def all_energies(self):
        energies = []
        for i in self.ionic_steps:
            en = float(i["energy"]["i"][1]["#text"])
            energies.append(en)
        return np.array(energies)

    @property
    def is_spin_polarized(self):
        if self.all_input_parameters["ISPIN"] == "2":
            return True
        else:
            return False

    @property
    def is_spin_orbit(self):
        if self.all_input_parameters["LSORBIT"] == "T":
            return True
        else:
            return False

    @property
    def all_structures(self):
        structs = []
        for i in self.ionic_steps:
            s = i["structure"]
            atoms = self.vrun_structure_to_atoms(s)
            structs.append(atoms)
        return structs

    @property
    def eigenvalues(self):
        nkpts = len(self.kpoints._kpoints)
        all_up_eigs = []
        all_dn_eigs = []
        if self.is_spin_polarized:
            for j in range(nkpts):
                eigs = np.array(
                    [
                        [float(jj) for jj in ii.split()]
                        for ii in (
                            self.ionic_steps[-1]["eigenvalues"]["array"]["set"]["set"][
                                0
                            ]
                        )["set"][j]["r"]
                    ]
                )
                all_up_eigs.append(eigs)
            for j in range(nkpts):
                eigs = np.array(
                    [
                        [float(jj) for jj in ii.split()]
                        for ii in (
                            self.ionic_steps[-1]["eigenvalues"]["array"]["set"]["set"][
                                1
                            ]
                        )["set"][j]["r"]
                    ]
                )
                all_dn_eigs.append(eigs)
        else:
            for j in range(nkpts):
                eigs = np.array(
                    [
                        [float(jj) for jj in ii.split()]
                        for ii in (
                            self.ionic_steps[-1]["eigenvalues"]["array"]["set"]["set"]
                        )["set"][j]["r"]
                    ]
                )
                all_up_eigs.append(eigs)
            all_dn_eigs = all_up_eigs

        all_up_eigs = np.array(all_up_eigs)
        all_dn_eigs = np.array(all_dn_eigs)
        return all_up_eigs, all_dn_eigs

    @property
    def all_forces(self):
        forces = []
        for m in self.ionic_steps:
            force = np.array(
                [[float(j) for j in i.split()] for i in m["varray"][0]["v"]]
            )

            forces.append(force)
        return np.array(forces)

    @property
    def all_stresses(self):
        stresses = []
        for m in self.ionic_steps:
            stress = np.array(
                [[float(j) for j in i.split()] for i in m["varray"][1]["v"]]
            )

            stresses.append(stress)
        return np.array(stresses)

    @property
    def all_input_parameters(self):
        d = OrderedDict()
        # import type
        for i in self._data["modeling"]["parameters"]["separator"]:
            for j, k in i.items():
                if j == "i":
                    for m in k:
                        if "#text" in m:
                            d[m["@name"]] = m["#text"]
                else:
                    if type(k) is list:
                        for n in k:
                            for p, q in n.items():
                                if p == "i":
                                    for r in q:
                                        if "#text" in r:
                                            d[r["@name"]] = r["#text"]
                                else:
                                    if type(q) is list:
                                        for s in q:
                                            if "#text" in s:
                                                d[s["@name"]] = s["#text"]
        return d

    @property
    def kpoints(self):
        kplist = np.array(
            [
                [float(j) for j in i.split()]
                for i in self._data["modeling"]["kpoints"]["varray"][0]["v"]
            ]
        )
        kpwt = np.array(
            [float(i) for i in self._data["modeling"]["kpoints"]["varray"][1]["v"]]
        )
        return Kpoints(kpoints=kplist, kpoints_weights=kpwt)

    def get_bandstructure(
        self, E_low=-4, E_high=4, spin=0, zero_efermi=True, kpoints_file_path="."
    ):
        try:
            f = open(kpoints_file_path, "r")
            lines = f.read().splitlines()
            f.close()
            kp_labels = []
            kp_labels_points = []
            for ii, i in enumerate(lines):
                if ii > 2:
                    tmp = i.split()
                    if len(tmp) == 5:
                        tmp = str("$") + str(tmp[4]) + str("$")
                        if len(kp_labels) == 0:
                            kp_labels.append(tmp)
                            kp_labels_points.append(ii - 3)
                        elif tmp != kp_labels[-1]:
                            kp_labels.append(tmp)
                            kp_labels_points.append(ii - 3)

        except:
            print("No K-points file found, still proceeding")
            pass

        tmp = 0.0
        if zero_efermi == True:
            tmp = float(self.efermi)
        plt.clf()
        for i, ii in enumerate(self.eigenvalues[spin][:, :, 0].T - tmp):
            plt.plot(ii, color="r")

        if self.is_spin_polarized:
            for i, ii in enumerate(self.eigenvalues[1][:, :, 0].T - tmp):
                plt.plot(ii, color="b")

        plt.ylim([E_low, E_high])
        plt.xticks(kp_labels_points, kp_labels)
        plt.xlim([0, len(self.kpoints._kpoints)])
        plt.xlabel(r"$\mathrm{Wave\ Vector}$")
        ylabel = (
            r"$\mathrm{E\ -\ E_f\ (eV)}$" if zero_efermi else r"$\mathrm{Energy\ (eV)}$"
        )
        plt.ylabel(ylabel)

        return plt

    @property
    def total_dos(self):
        energies = []
        spin_up = []
        spin_up_data = np.array(
            [
                [float(j) for j in i.split()]
                for i in self.ionic_steps[-1]["dos"]["total"]["array"]["set"]["set"][0][
                    "r"
                ]
            ]
        )
        energies = spin_up_data[:, 0]
        spin_up = spin_up_data[:, 1]
        if self.is_spin_polarized:
            spin_dn = []
            spin_dn_data = np.array(
                [
                    [float(j) for j in i.split()]
                    for i in self.ionic_steps[-1]["dos"]["total"]["array"]["set"][
                        "set"
                    ][1]["r"]
                ]
            )
            spin_dn = -1 * spin_dn_data[:, 1]
        return energies, spin_up, spin_dn


class Oszicar(object):
    def __init__(self, filename, data={}):
        """
        Class handling VASP OSZICAR file data
        """
        self.filename = filename
        self.data = data
        if self.data == {}:
            f = open(filename, "r")
            lines = f.read().splitlines()
            f.close()
            self.data = lines

    @property
    def magnetic_moment(self):
        return self.ionic_steps[-1][-1]

    @property
    def ionic_steps(self):
        ionic_data = []
        for i in self.data:
            if "E0" in i:
                ionic_data.append(i.split())
        return ionic_data

    @property
    def electronic_steps(self):
        electronic_data = []
        for i in self.data:
            if "E0" not in i and "d eps" not in i:
                electronic_data.append(i.split())
        return electronic_data


class Outcar(object):
    def __init__(self, filename, data={}):
        """
        Class handling VASP OUTCAR file data
        """
        self.filename = filename
        self.data = data
        if self.data == {}:
            f = open(filename, "r")
            lines = f.read().splitlines()
            f.close()
            self.data = lines

    @property
    def nions(self):
        for i in self.data:
            if "NIONS =" in i:
                n_ions = int(i.split()[-1])
                return n_ions

    @property
    def phonon_eigenvalues(self):
        # Thz values
        lines = self.data
        vals = []
        for i, ii in enumerate(lines):
            if "meV" in ii:
                tmp = float(ii.split()[-8])
                if "f/i" in ii:
                    tmp = tmp * -1
                vals.append(tmp)
        return np.array(vals, dtype="float")

    @property
    def converged(self):

        cnvg = False
        try:
            lines = self.data
            cnvg = False
            for i in lines:
                if "General timing and accounting informations for this job" in i:
                    cnvg = True
            # print fil,cnvg
        except:
            pass
        return cnvg

    @property
    def efg_tensor_diag(self):
        nions = self.nions
        for ii, i in enumerate(self.data):
            if "Electric field gradients after diagonalization" in i:
                tmp = ii
        arr = self.data[tmp + 5 : tmp + 5 + nions]
        efg_arr = []
        for i in arr:
            tmp = [i.split()[1], i.split()[2], i.split()[3], i.split()[4]]
            efg_arr.append(tmp)
        efg_arr = np.array(efg_arr, dtype="float")
        return efg_arr

    @property
    def quad_mom(self):
        nions = self.nions
        for ii, i in enumerate(self.data):
            if "Q  : nuclear electric quadrupole moment in mb (millibarn)" in i:
                tmp = ii
        arr = self.data[tmp + 4 : tmp + 4 + nions]
        quad_arr = []
        for i in arr:
            tmp = [i.split()[1], i.split()[2], i.split()[3]]
            quad_arr.append(tmp)
        quad_arr = np.array(quad_arr, dtype="float")
        return quad_arr

    @property
    def piezoelectric_tensor(self):
        lines = self.data
        ionic_piezo = []
        total_piezo = []
        for ii, i in enumerate(lines):
            if "PIEZOELECTRIC TENSOR" in i and "(C/m^2)" in i and "field" in i:
                if "IONIC" in i:
                    ionic_piezo.append(lines[ii + 3].split()[1:7])
                    ionic_piezo.append(lines[ii + 4].split()[1:7])
                    ionic_piezo.append(lines[ii + 5].split()[1:7])
                else:
                    total_piezo.append(lines[ii + 3].split()[1:7])
                    total_piezo.append(lines[ii + 4].split()[1:7])
                    total_piezo.append(lines[ii + 5].split()[1:7])
        ionic_piezo = np.array(ionic_piezo, dtype="float")
        total_piezo = np.array(total_piezo, dtype="float")
        return ionic_piezo, total_piezo

    def elastic_props(self, atoms=None, vacuum=False):
        """
            Obtain elastic tensor and calculate related properties
            Args:
                outcar: OUTCAR file path
                vacuum: whether the structure has vaccum such as 2D materials
                for vacuum structures bulk and shear mod. needs extra attenstion
                and elastic tensor are in Nm^-1 rather than GPa
            Returns:
                  info: data for elastic tensor (in string and object representation), bulk, shear modulus, and phonon modes
            """
        ratio_c = 1.0
        if vacuum == True:
            ratio_c = 0.1 * float(
                abs(atoms.lattice_mat[2][2])
            )  # *(10**9)*(10**-10) #N/m unit
        KV = "na"
        GV = "na"
        spin = "na"
        info = {}
        v = open(self.filename, "r")
        lines = v.read().splitlines()
        c = np.empty((6, 6), dtype=float)
        for i, line in enumerate(lines):
            if "TOTAL ELASTIC MODULI (kBar)" in line:
                c11 = lines[i + 3].split()[1]
                c12 = lines[i + 3].split()[2]
                c13 = lines[i + 3].split()[3]
                c14 = lines[i + 3].split()[4]
                c15 = lines[i + 3].split()[5]
                c16 = lines[i + 3].split()[6]
                c21 = lines[i + 4].split()[1]
                c22 = lines[i + 4].split()[2]
                c23 = lines[i + 4].split()[3]
                c24 = lines[i + 4].split()[4]
                c25 = lines[i + 4].split()[5]
                c26 = lines[i + 4].split()[6]
                c31 = lines[i + 5].split()[1]
                c32 = lines[i + 5].split()[2]
                c33 = lines[i + 5].split()[3]
                c34 = lines[i + 5].split()[4]
                c35 = lines[i + 5].split()[5]
                c36 = lines[i + 5].split()[6]
                c41 = lines[i + 6].split()[1]
                c42 = lines[i + 6].split()[2]
                c43 = lines[i + 6].split()[3]
                c44 = lines[i + 6].split()[4]
                c45 = lines[i + 6].split()[5]
                c46 = lines[i + 6].split()[6]
                c51 = lines[i + 7].split()[1]
                c52 = lines[i + 7].split()[2]
                c53 = lines[i + 7].split()[3]
                c54 = lines[i + 7].split()[4]
                c55 = lines[i + 7].split()[5]
                c56 = lines[i + 7].split()[6]
                c61 = lines[i + 8].split()[1]
                c62 = lines[i + 8].split()[2]
                c63 = lines[i + 8].split()[3]
                c64 = lines[i + 8].split()[4]
                c65 = lines[i + 8].split()[5]
                c66 = lines[i + 8].split()[6]
                c[0][0] = round(ratio_c * float(c11) / float(10), 1)
                c[0][1] = round(ratio_c * float(c12) / float(10), 1)
                c[0][2] = round(ratio_c * float(c13) / float(10), 1)
                c[0][3] = round(ratio_c * float(c14) / float(10), 1)
                c[0][4] = round(ratio_c * float(c15) / float(10), 1)
                c[0][5] = round(ratio_c * float(c16) / float(10), 1)
                c[1][0] = round(ratio_c * float(c21) / float(10), 1)
                c[1][1] = round(ratio_c * float(c22) / float(10), 1)
                c[1][2] = round(ratio_c * float(c23) / float(10), 1)
                c[1][3] = round(ratio_c * float(c24) / float(10), 1)
                c[1][4] = round(ratio_c * float(c25) / float(10), 1)
                c[1][5] = round(ratio_c * float(c26) / float(10), 1)
                c[2][0] = round(float(c31) / float(10), 1)
                c[2][1] = round(float(c32) / float(10), 1)
                c[2][2] = round(float(c33) / float(10), 1)
                c[2][3] = round(float(c34) / float(10), 1)
                c[2][4] = round(float(c35) / float(10), 1)
                c[2][5] = round(float(c36) / float(10), 1)
                c[3][0] = round(float(c41) / float(10), 1)
                c[3][1] = round(float(c42) / float(10), 1)
                c[3][2] = round(float(c43) / float(10), 1)
                c[3][3] = round(float(c44) / float(10), 1)
                c[3][4] = round(float(c45) / float(10), 1)
                c[3][5] = round(float(c46) / float(10), 1)
                c[4][0] = round(float(c51) / float(10), 1)
                c[4][1] = round(float(c52) / float(10), 1)
                c[4][2] = round(float(c53) / float(10), 1)
                c[4][3] = round(float(c54) / float(10), 1)
                c[4][4] = round(float(c55) / float(10), 1)
                c[4][5] = round(float(c56) / float(10), 1)
                c[5][0] = round(float(c61) / float(10), 1)
                c[5][1] = round(float(c62) / float(10), 1)
                c[5][2] = round(float(c63) / float(10), 1)
                c[5][3] = round(float(c64) / float(10), 1)
                c[5][4] = round(float(c65) / float(10), 1)
                c[5][5] = round(float(c66) / float(10), 1)
                KV = float(
                    (c[0][0] + c[1][1] + c[2][2]) + 2 * (c[0][1] + c[1][2] + c[2][0])
                ) / float(9)
                GV = float(
                    (c[0][0] + c[1][1] + c[2][2])
                    - (c[0][1] + c[1][2] + c[2][0])
                    + 3 * (c[3][3] + c[4][4] + c[5][5])
                ) / float(15)
                KV = round(KV, 3)
                GV = round(GV, 3)
                break
        v.close()

        modes = []
        try:
            for i in lines:
                if "cm-1" in i and "meV" in i:

                    mod = float(i.split()[-4])
                    if "f/i" in i:
                        mod = mod * -1
                    if mod not in modes:
                        modes.append(float(mod))
        except:
            pass

        info["cij"] = c.tolist()
        info["KV"] = KV
        info["GV"] = GV
        info["modes"] = modes

        return info


class Waveder(object):
    """
    Class for reading a WAVEDER file.
    The LOPTICS tag produces a WAVEDER file.
    The WAVEDER contains the derivative of the orbitals with respect to k.
    """

    def __init__(self, filename, gamma_only=False):
        """
        Args:
            filename: Name of file containing WAVEDER.
        """
        with open(filename, "rb") as fp:

            def readData(dtype):
                """ Read records from Fortran binary file and convert to
                np.array of given dtype. """
                data = b""
                while 1:
                    prefix = np.fromfile(fp, dtype=np.int32, count=1)[0]
                    data += fp.read(abs(prefix))
                    suffix = np.fromfile(fp, dtype=np.int32, count=1)[0]
                    if abs(prefix) - abs(suffix):
                        raise RuntimeError(
                            "Read wrong amount of bytes.\n"
                            "Expected: %d, read: %d, suffix: %d."
                            % (prefix, len(data), suffix)
                        )
                    if prefix > 0:
                        break
                return np.frombuffer(data, dtype=dtype)

            nbands, nelect, nk, ispin = readData(np.int32)
            _ = readData(np.float)  # nodes_in_dielectric_function
            _ = readData(np.float)  # wplasmon
            if gamma_only:
                cder = readData(np.float)
            else:
                cder = readData(np.complex64)

            cder_data = cder.reshape((3, ispin, nk, nelect, nbands)).T

            self._cder_data = cder_data
            self._nkpoints = nk
            self._ispin = ispin
            self._nelect = nelect
            self._nbands = nbands

    @property
    def cder_data(self):
        """
        Returns the orbital derivative between states
        """
        return self._cder_data

    @property
    def nbands(self):
        """
        Returns the number of bands in the calculation
        """
        return self._nbands

    @property
    def nkpoints(self):
        """
        Returns the number of k-points in the calculation
        """
        return self._nkpoints

    @property
    def nelect(self):
        """
        Returns the number of electrons in the calculation
        """
        return self._nelect

    def get_orbital_derivative_between_states(
        self, band_i, band_j, kpoint, spin, cart_dir
    ):
        """
        Method returning a value
        between bands band_i and band_j for k-point index, spin-channel and cartesian direction.
        Args:
            band_i (Integer): Index of band i
            band_j (Integer): Index of band j
            kpoint (Integer): Index of k-point
            spin   (Integer): Index of spin-channel (0 or 1)
            cart_dir (Integer): Index of cartesian direction (0,1,2)
        Returns:
            a float value
        """
        if (
            band_i < 0
            or band_i > self.nbands - 1
            or band_j < 0
            or band_j > self.nelect - 1
        ):
            raise ValueError("Band index out of bounds")
        if kpoint > self.nkpoints:
            raise ValueError("K-point index out of bounds")
        if cart_dir > 2 or cart_dir < 0:
            raise ValueError("cart_dir index out of bounds")

        return self._cder_data[band_i, band_j, kpoint, spin, cart_dir]


class Wavecar(object):
    """
    Class for VASP Pseudowavefunction stored in WAVECAR

    The format of VASP WAVECAR, as shown in
        http://www.andrew.cmu.edu/user/feenstra/wavetrans/
    is:
        Record-length #spin components RTAG(a value specifying the precision)
        #k-points #bands ENCUT(maximum energy for plane waves)
        LatVec-A
        LatVec-B
        LatVec-C
        Loop over spin
           Loop over k-points
              #plane waves, k vector
              Loop over bands
                 band energy, band occupation
              End loop over bands
              Loop over bands
                 Loop over plane waves
                    Plane-wave coefficient
                 End loop over plane waves
              End loop over bands
           End loop over k-points
        End loop over spin
    """

    def __init__(
        self,
        filename="WAVECAR",
        lsorbit=False,
        lgamma=False,
        recl=None,
        nspin=None,
        rtag=None,
        nkpts=None,
        nbands=None,
        encut=None,
        lattice_mat=None,
        nplws=None,
        wfc=None,
        efermi=None,
        kvecs=None,
        energies=None,
        occs=None,
        gvec=None,
    ):
        """
        Initialization.
        """

        self._filename = filename
        self._lsoc = lsorbit
        self._lgam = lgamma
        self._recl = recl
        self._nspin = nspin
        self._rtag = rtag
        self._nkpts = nkpts
        self._nbands = nbands
        self._encut = encut
        self._lattice_mat = lattice_mat
        self._nplws = nplws
        self._efermi = efermi
        self._kvecs = kvecs
        self._energies = energies
        self._occs = occs
        self._gvec = gvec
        self._lattice_mat = lattice_mat

        assert not (lsorbit and lgamma), "The two settings conflict!"

        try:
            self._wfc = open(self._filename, "rb")
        except:
            raise IOError("Failed to open %s" % self._fname)

        # read the basic information
        self.readWFHeader()
        # read the band information
        self.readWFBand()

        if self._lsoc:
            assert self._nspin == 1, "NSPIN = 1 for noncollinear version WAVECAR!"

    def isSocWfc(self):
        """
        Is the WAVECAR from an SOC calculation?
        """
        return True if self._lsoc else False

    def readWFHeader(self):
        """
        Read the system information from WAVECAR, which is written in the first
        two record.

        rec1: recl, nspin, rtag
        rec2: nkpts, nbands, encut, ((cell(i,j) i=1, 3), j=1, 3)
        """

        # goto the start of the file and read the first record
        self._wfc.seek(0)
        self._recl, self._nspin, self._rtag = np.array(
            np.fromfile(self._wfc, dtype=np.float, count=3), dtype=int
        )
        self._WFPrec = self.setWFPrec()
        # the second record
        self._wfc.seek(self._recl)
        dump = np.fromfile(self._wfc, dtype=np.float, count=12)

        self._nkpts = int(dump[0])  # No. of k-points
        self._nbands = int(dump[1])  # No. of bands
        self._encut = dump[2]  # Energy cutoff
        self._lattice_mat = dump[3:].reshape((3, 3))  # real space supercell basis
        self._Omega = np.linalg.det(self._lattice_mat)  # real space supercell volume
        self._Bcell = np.linalg.inv(
            self._lattice_mat
        ).T  # reciprocal space supercell volume

        # Minimum FFT grid size
        Anorm = np.linalg.norm(self._lattice_mat, axis=1)
        CUTOF = np.ceil(np.sqrt(self._encut / RYTOEV) / (TPI / (Anorm / AUTOA)))
        self._ngrid = np.array(2 * CUTOF + 1, dtype=int)

    def setWFPrec(self):
        """
        Set wavefunction coefficients precision:
            TAG = 45200: single precision complex, np.complex64, or complex(qs)
            TAG = 45210: double precision complex, np.complex128, or complex(q)
        """
        if self._rtag == 45200:
            return np.complex64
        elif self._rtag == 45210:
            return np.complex128
        elif self._rtag == 53300:
            raise ValueError("VASP5 WAVECAR format, not implemented yet")
        elif self._rtag == 53310:
            raise ValueError(
                "VASP5 WAVECAR format with double precision "
                + "coefficients, not implemented yet"
            )
        else:
            raise ValueError("Invalid TAG values: {}".format(self._rtag))

    def readWFBand(self, ispin=1, ikpt=1, iband=1):
        """
        Extract KS energies and Fermi occupations from WAVECAR.
        """

        self._nplws = np.zeros(self._nkpts, dtype=int)
        self._kvecs = np.zeros((self._nkpts, 3), dtype=float)
        self._energies = np.zeros((self._nspin, self._nkpts, self._nbands), dtype=float)
        self._occs = np.zeros((self._nspin, self._nkpts, self._nbands), dtype=float)

        for ii in range(self._nspin):
            for jj in range(self._nkpts):
                rec = self.whereRec(ii + 1, jj + 1, 1) - 1
                self._wfc.seek(rec * self._recl)
                dump = np.fromfile(
                    self._wfc, dtype=np.float, count=4 + 3 * self._nbands
                )
                if ii == 0:
                    self._nplws[jj] = int(dump[0])
                    self._kvecs[jj] = dump[1:4]
                dump = dump[4:].reshape((-1, 3))
                self._energies[ii, jj, :] = dump[:, 0]
                self._occs[ii, jj, :] = dump[:, 2]

        if self._nkpts > 1:
            tmp = np.linalg.norm(
                np.dot(np.diff(self._kvecs, axis=0), self._Bcell), axis=1
            )
            self._kpath = np.concatenate(([0], np.cumsum(tmp)))
        else:
            self._kpath = None
        return self._kpath, self._energies

    def gvectors(self, ikpt=1):
        """
        Generate the G-vectors that satisfies the following relation
            (G + k)**2 / 2 < ENCUT
        """
        assert 1 <= ikpt <= self._nkpts, "Invalid kpoint index!"

        kvec = self._kvecs[ikpt - 1]
        # fx, fy, fz = [fftfreq(n) * n for n in self._ngrid]
        # fftfreq in scipy.fftpack is a little different with VASP frequencies
        fx = [
            ii if ii < self._ngrid[0] / 2 + 1 else ii - self._ngrid[0]
            for ii in range(self._ngrid[0])
        ]
        fy = [
            jj if jj < self._ngrid[1] / 2 + 1 else jj - self._ngrid[1]
            for jj in range(self._ngrid[1])
        ]
        fz = [
            kk if kk < self._ngrid[2] / 2 + 1 else kk - self._ngrid[2]
            for kk in range(self._ngrid[2])
        ]
        if self._lgam:
            # parallel gamma version of VASP WAVECAR exclude some planewave
            # components, -DwNGZHalf
            kgrid = np.array(
                [
                    (fx[ii], fy[jj], fz[kk])
                    for kk in range(self._ngrid[2])
                    for jj in range(self._ngrid[1])
                    for ii in range(self._ngrid[0])
                    if (
                        (fz[kk] > 0)
                        or (fz[kk] == 0 and fy[jj] > 0)
                        or (fz[kk] == 0 and fy[jj] == 0 and fx[ii] >= 0)
                    )
                ],
                dtype=float,
            )
        else:
            kgrid = np.array(
                [
                    (fx[ii], fy[jj], fz[kk])
                    for kk in range(self._ngrid[2])
                    for jj in range(self._ngrid[1])
                    for ii in range(self._ngrid[0])
                ],
                dtype=float,
            )

        # Kinetic_Energy = (G + k)**2 / 2
        # HSQDTM    =  hbar**2/(2*ELECTRON MASS)
        KENERGY = (
            HSQDTM
            * np.linalg.norm(
                np.dot(kgrid + kvec[np.newaxis, :], TPI * self._Bcell), axis=1
            )
            ** 2
        )
        # find Gvectors where (G + k)**2 / 2 < ENCUT
        Gvec = kgrid[np.where(KENERGY < self._encut)[0]]

        if self._lsoc:
            assert Gvec.shape[0] == self._nplws[ikpt - 1] / 2, (
                "No. of planewaves not consistent for an SOC WAVECAR! %d %d %d"
                % (Gvec.shape[0], self._nplws[ikpt - 1], np.prod(self._ngrid))
            )
        else:
            assert Gvec.shape[0] == self._nplws[ikpt - 1], (
                "No. of planewaves not consistent! %d %d %d"
                % (Gvec.shape[0], self._nplws[ikpt - 1], np.prod(self._ngrid))
            )
        self._gvec = np.asarray(Gvec, dtype=int)

        return np.asarray(Gvec, dtype=int)

    def readBandCoeff(self, ispin=1, ikpt=1, iband=1, norm=False):
        """
        Read the planewave coefficients of specified KS states.
        """

        self.checkIndex(ispin, ikpt, iband)

        rec = self.whereRec(ispin, ikpt, iband)
        self._wfc.seek(rec * self._recl)

        nplw = self._nplws[ikpt - 1]
        dump = np.fromfile(self._wfc, dtype=self._WFPrec, count=nplw)

        cg = np.asarray(dump, dtype=np.complex128)
        if norm:
            cg /= np.linalg.norm(cg)
        return cg

    def whereRec(self, ispin=1, ikpt=1, iband=1):
        """
        Return the rec position for specified KS state.
        """

        self.checkIndex(ispin, ikpt, iband)

        rec = (
            2
            + (ispin - 1) * self._nkpts * (self._nbands + 1)
            + (ikpt - 1) * (self._nbands + 1)
            + iband
        )
        return rec

    def checkIndex(self, ispin, ikpt, iband):
        """
        Check if the index is valid!
        """
        assert 1 <= ispin <= self._nspin, "Invalid spin index!"
        assert 1 <= ikpt <= self._nkpts, "Invalid kpoint index!"
        assert 1 <= iband <= self._nbands, "Invalid band index!"


"""
if __name__ == "__main__":
    
    o=Outcar('../../tests/testfiles/io/vasp/OUTCAR.JVASP-39')
    print (round(o.phonon_eigenvalues[2],2))
    #print (o.born_eff_charg)
    sys.exit()
    v=Vasprun('../../tests/testfiles/io/vasp/vasprun.xml.JVASP-39')
    dfpt=round(v.dfpt_data['born_charges'][0][0][0],2)
    print (dfpt)
    o = Outcar(
        "/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-134_bulk_PBEBO/MAIN-ELASTIC-bulk@mp-134/OUTCAR"
    )
    o = Outcar("/rk2/knc6/EFG/JVASP-1429/MAIN-RELAX-bulk@mp_2964/OUTCAR")
    o = Outcar(
        "/rk2/knc6/EFG/JVASP-4798_mp-12597_PBEBO/MAIN-LEFG-JVASP-4798_mp-12597/OUTCAR"
    )
    print("NIONS2", o.efg_tensor_diag)
    print("NIONS2", o.efg_tensor_diag[:, 0:3])
    # print ('NIONS2',o.quad_mom)
    o = Outcar(
        "/rk2/knc6/JARVIS-DFT/Solar-Semi/mp-5986_bulk_PBEBO/MAIN-ELASTIC-bulk@mp_5986/OUTCAR"
    )
    filename = "/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/CHGCAR"
    x = Chgcar(filename).chg
    print(len(x), x[0])
    # print (o.elastic_props())
    import sys

    sys.exit()
    v = Vasprun(
        "/cluster/users/knc6/justback/Interfa/DFT/JVASP-649_JVASP-76195_PBEBO/RELAXSOCPBEBAND/vasprun.xml"
    )
    # v=Vasprun('/cluster/users/knc6/justback/Interfa/DFT/JVASP-649_JVASP-76195_PBEBO/RELAXPBEBAND/vasprun.xml')

    p = v.get_bandstructure(
        kpoints_file_path="/cluster/users/knc6/justback/Interfa/DFT/JVASP-649_JVASP-76195_PBEBO/RELAXPBEBAND/KPOINTS"
    )
    p.savefig("PBEBAND2.png")

    sys.exit()
    oz = Oszicar(
        "/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/OSZICAR"
    )
    print(oz.electronic_steps)
    # c=Chgcar()
    # filename = "/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/CHGCAR"
    # filename = "/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/LOCPOT"
    # x = Chgcar(filename).chg
    # print(len(x), x[0])
    # y = VaspChargeDensity(filename).chg
    # print(len(y), y[0])
    # if x[0].all() == y[0].all():
    #    print("JACKPOT")
    # filename='/rk2/knc6/Chern3D/JVASP-1067_mp-541837_PBEBO/MAIN-SOC-bulk@JVASP-1067_mp-541837/CHGCAR'
    # x=Chgcar(filename)
    # filename='/rk2/knc6/Chern3D/JVASP-1067_mp-541837_PBEBO/MAIN-MAGSCF-bulk@JVASP-1067_mp-541837/CHGCAR'
    # x=Chgcar(filename)
    # filename='/rk2/knc6/Chern/JVASP-13600_mp-28208_PBEBO/MAIN-MAGSCF-bulk@JVASP-13600_mp-28208/CHGCAR'
    # x=Chgcar(filename)
    wf_so = "/rk2/knc6/Chern3D/JVASP-1067_mp-541837_PBEBO/MAIN-SOCSCFBAND-bulk@JVASP-1067_mp-541837/WAVECAR"
    so = Wavecar(filename=wf_so, lsorbit=True)
    print("Read SO")
    sys.exit()
    # filename = "/rk2/knc6/Chern3D/JVASP-1067_mp-541837_PBEBO/MAIN-SOCSCFBAND-bulk@JVASP-1067_mp-541837/vasprun.xml"
    # filename='/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-134_bulk_PBEBO/MAIN-RELAX-bulk@mp-134/vasprun.xml'
    # filename = '/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/vasprun.xml'
    filename = "/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-OPTICS-bulk@mp-149/vasprun.xml"
    filename = "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-RELAX-bulk@mp_541837/vasprun.xml"
    filename = "/rk2/knc6/JARVIS-DFT/Bulk-less5mp1/mp-1077840_PBEBO/MAIN-RELAX-bulk@mp_1077840/vasprun.xml"  # ????
    v = Vasprun(filename=filename)
    # print (v.is_spin_orbit,v.is_spin_polarized)

    print("dir", v.get_dir_gap, "indir", v.get_indir_gap)
    # print ('reals',v.avg_absorption_coefficient)
    from pymatgen.io.vasp.outputs import Vasprun as pmg

    pmgv = pmg(filename)
    bandstructure = pmgv.get_band_structure()
    print("pmggap", bandstructure.get_direct_band_gap(), bandstructure.get_band_gap())
    import sys

    filename = "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-RELAX-bulk@mp_541837/OSZICAR"
    osz = Oszicar(filename=filename)
    # print (osz.ionic_steps(),osz.magnetic_moment())
    import sys

    # print (v._filename, v.final_energy)
    # print (v._filename,'elements', v.elements)
    # print(v._filename, "structs", v.all_structures)
    # print ('pmg',Vasprun(v._filename).final_energy)
    filename = "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-BAND-bulk@mp_541837/vasprun.xml"
    v = Vasprun(filename=filename)
    # print(v._filename, "structs", v.all_structures)
    # print (v._filename,v.final_energy)
    # print (v._filename,v.elements)
    # print ('pmg',Vasprun(v._filename).final_energy)
    # print (v._data.keys())
    # print ()
    # print ()
    # print (v._data['modeling'].keys())
    # print ()
    # print ()
    ##'generator', 'incar', 'kpoints', 'parameters', 'atominfo', 'structure', 'calculation'
    # print ('incar',v._data['modeling']['incar'])
    # print ()
    # print ()
    # print ('kpoints',v._data['modeling']['kpoints'])
    # print ()
    # print ()
    # print ('atominfo',v._data['modeling']['atominfo'])
    # print ()
    # print ()
    # print ('parameters',v._data['modeling']['parameters'])
    # print ()
    # print ()
    # print ('structure',v._data['modeling']['structure'])
    # print ()
    # print ()
    # print ('calculation',v._data['modeling']['calculation'].keys())
    # print ()
    # print ()
    ##calculation odict_keys(['scstep', 'structure', 'varray', 'energy', 'time', 'eigenvalues', 'separator', 'dos', 'projected'])
"""
