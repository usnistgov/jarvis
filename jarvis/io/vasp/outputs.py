"""Modulet for analzing VASP outputs."""

from scipy.constants import physical_constants
from scipy.constants import speed_of_light
from jarvis.core.atoms import Atoms
import numpy as np
from collections import OrderedDict
from jarvis.core.specie import Specie
import xmltodict
from jarvis.core.kpoints import Kpoints3D as Kpoints
from jarvis.io.vasp.inputs import Poscar
from matplotlib import pyplot as plt
from jarvis.core.utils import rec_dict
from jarvis.core.utils import recast_array_on_uniq_array_elements
import scipy.signal as ss

RYTOEV = 13.605826
AUTOA = 0.529177249
TPI = 2 * np.pi
HSQDTM = RYTOEV * AUTOA * AUTOA
plt.switch_backend("agg")


class Chgcar(object):
    """Class handling VASP CHGCAR file data."""

    def __init__(
        self,
        filename="",
        atoms=None,
        chg=[],
        chgdif=None,
        aug=None,
        augdiff=None,
        dim=None,
        nsets=1,
    ):
        """
        Contain CHGCAR data.

        Requires following arguments, but some are optional.

        Args:
           filename: path of CHGCAR

           atoms: Atoms

           chg: charge density

           chgdif: difference in up and down spin

           aug: PAW augmentation

           augdiff: augmentation difference

           nsets: number of CHG sets.
        """
        self.filename = filename
        self.atoms = atoms
        self.chg = chg
        self.dim = dim
        self.chgdif = chgdif
        self.aug = aug
        self.augdiff = augdiff
        self.nsets = nsets
        if self.atoms is None:
            self.read_file()

    def to_dict(self):
        """Convert to a dictionary."""
        d = OrderedDict()
        d["filename"] = self.filename
        d["atoms"] = self.atoms.to_dict()
        d["chg"] = self.chg
        d["dim"] = self.dim
        d["chgdif"] = self.chgdif
        d["aug"] = self.aug
        d["augdiff"] = self.augdiff
        d["nsets"] = self.nsets
        if self.atoms is not None:
            d["atoms"] = self.atoms.to_dict()
        else:
            d["atoms"] = self.atoms
        return d

    @classmethod
    def from_dict(self, d={}):
        """Construct class from a dictionary."""
        if d["atoms"] is not None:
            atoms = Atoms.from_dict(d["atoms"])
        else:
            atoms = None
        return Chgcar(
            filename=d["filename"],
            atoms=atoms,
            chg=d["chg"],
            dim=d["dim"],
            chgdif=d["chgdif"],
            aug=d["aug"],
            augdiff=d["augdiff"],
            nsets=d["nsets"],
        )

    def is_spin_polarized(self):
        """Check if the calculations is spin-polarized, ISPIN=2."""
        if self.nsets == 2:
            return True
        else:
            return False

    def is_spin_orbit(self):
        """Check if the calculations is spin-orbit, LSORBIT=T."""
        if self.nsets == 4:
            return True
        else:
            return False

    def read_file(self):
        """Read CHGCAR."""
        f = open(self.filename, "r")
        lines = f.read()
        f.close()
        self.atoms = Poscar.from_string(lines).atoms
        volume = self.atoms.volume
        text = lines.splitlines()
        ng_line = text[self.atoms.num_atoms + 9]
        ng = [int(j) for j in ng_line.split()]
        self.dim = np.array(ng)
        found = self.atoms.num_atoms + 8
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
        self.chg = np.array(self.chg)

    def chg_set(self, text, start, end, volume, ng):
        """Return CHGCAR sets."""
        lines_0 = text[start:end]
        tmp = []
        for i in lines_0:
            for j in i.split():
                if j != "":
                    tmp.append(float(j))
        tmp = np.array(tmp).reshape(ng)
        tmp = tmp / volume
        return tmp


class Locpot(Chgcar):
    """Read LOCPOT files."""

    def vac_potential(
        self, direction="X", Ef=0, filename="Avg.png", plot=True
    ):
        """Calculate vacuum potential used in work-function calculation."""
        atoms = self.atoms
        cell = atoms.lattice_mat
        chg = (self.chg[-1].T) * atoms.volume
        latticelength = np.dot(cell, cell.T).diagonal()
        latticelength = latticelength ** 0.5
        ngridpts = np.array(chg.shape)
        # totgridpts = ngridpts.prod()

        if direction == "X":
            idir = 0
            a = 1
            b = 2
        elif direction == "Y":
            a = 0
            idir = 1
            b = 2
        else:
            a = 0
            b = 1
            idir = 2
        a = (idir + 1) % 3
        b = (idir + 2) % 3
        average = np.zeros(ngridpts[idir], np.float)
        for ipt in range(ngridpts[idir]):
            if direction == "X":
                average[ipt] = chg[ipt, :, :].sum()
            elif direction == "Y":
                average[ipt] = chg[:, ipt, :].sum()
            else:
                average[ipt] = chg[:, :, ipt].sum()
        average /= ngridpts[a] * ngridpts[b]
        xdiff = latticelength[idir] / float(ngridpts[idir] - 1)
        xs = []
        ys = []
        for i in range(ngridpts[idir]):
            x = i * xdiff
            xs.append(x)
            ys.append(average[i])

        avg_max = max(average)

        dif = float(avg_max) - float(Ef)
        if plt:
            plt.xlabel("z (Angstrom)")
            plt.plot(xs, ys, "-", linewidth=2, markersize=10)
            horiz_line_data = np.array([avg_max for i in range(len(xs))])
            plt.plot(xs, horiz_line_data, "-")
            horiz_line_data = np.array([Ef for i in range(len(xs))])
            plt.plot(xs, horiz_line_data, "-")
            plt.ylabel("Potential (eV)")
            ax = plt.gca()
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            plt.title(
                str("Energy difference ")
                + str(round(float(dif), 3))
                + str(" eV"),
                fontsize=26,
            )
            plt.tight_layout()

            plt.savefig(filename)
            plt.close()

        print("Ef,max,wf=", Ef, avg_max, dif)
        return avg_max, dif


class Oszicar(object):
    """Construct Oszicar object."""

    def __init__(self, filename, data={}):
        """Initialize with filename."""
        self.filename = filename
        self.data = data
        if self.data == {}:
            f = open(filename, "r")
            lines = f.read().splitlines()
            f.close()
            self.data = lines

    @classmethod
    def from_dict(self, d={}):
        """Construct class from a dictionary."""
        return Oszicar(filename=d["filename"], data=d["data"])

    def to_dict(self):
        """Convert class to a dictionary."""
        d = OrderedDict()
        d["filename"] = self.filename
        d["data"] = self.data
        return d

    @property
    def magnetic_moment(self):
        """Get magnetic moment."""
        # return self.ionic_steps[-1][-1]
        match = -1
        for i, ii in enumerate(self.ionic_steps[-1]):
            if ii == "mag=":
                match = i + 1
        return np.array([float(j) for j in self.ionic_steps[-1][match:]])

    @property
    def ionic_steps(self):
        """Get all ionic steps realted data."""
        ionic_data = []
        for i in self.data:
            if "E0" in i:
                ionic_data.append(i.split())
        return ionic_data

    @property
    def electronic_steps(self):
        """Get all electronic steps realted data."""
        electronic_data = []
        for i in self.data:
            if "E0" not in i and "d eps" not in i:
                electronic_data.append(i.split())
        return electronic_data


class Outcar(object):
    """Construct OUTCAR object."""

    def __init__(self, filename, data={}):
        """Intialize with filename."""
        self.filename = filename
        self.data = data
        if self.data == {}:
            f = open(filename, "r")
            lines = f.read().splitlines()
            f.close()
            self.data = lines

    @classmethod
    def from_dict(self, d={}):
        """Construct class from a dictionary."""
        return Outcar(filename=d["filename"], data=d["data"])

    def to_dict(self):
        """Convert class to a dictionary."""
        d = OrderedDict()
        d["filename"] = self.filename
        d["data"] = self.data
        return d

    @property
    def nions(self):
        """Get number of ions."""
        for i in self.data:
            if "NIONS =" in i:
                n_ions = int(i.split()[-1])
                return n_ions

    @property
    def nbands(self):
        """Get number of bands."""
        for i in self.data:
            if "NBANDS=" in i:
                nbands = int(i.split()[-1])
                return nbands

    @property
    def nelect(self):
        """Get number of electrons."""
        for i in self.data:
            if "NELECT" in i:
                nelect = int(float(i.split()[2]))
                return nelect

    @property
    def phonon_eigenvalues(self):
        """Get phonon eigenvalues."""
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
        """Check if calculation is converged."""
        cnvg = False
        try:
            lines = self.data
            cnvg = False
            for i in lines:
                if (
                    "General timing and accounting informations for this job"
                    in i
                ):
                    cnvg = True
                if "VASP will stop now." in i:
                    cnvg = True
        except Exception:
            pass
        return cnvg

    #  @property
    def efg_tensor_diag(self, std_conv=True, prec=3):
        """Get diagonalized electric field gradient tensor."""
        # std_conv: |Vzz|>=|Vyy|>=|Vxx|, eta=(Vxx-Vyy)/Vzz
        # Note: VASP uses: |Vzz|>=|Vxx|>=|Vyy|, eta=(Vyy-Vxx)/Vzz
        # quadrupolar parameter, Cq=e*Q*V_zz/h
        nions = self.nions
        for ii, i in enumerate(self.data):
            if "Electric field gradients after diagonalization" in i:
                tmp = ii
        arr = self.data[tmp + 5 : tmp + 5 + nions]
        efg_arr = []
        for i in arr:
            if std_conv:
                Vzz = round(float(i.split()[3]), prec)
                Vyy = round(float(i.split()[1]), prec)
                Vxx = round(float(i.split()[2]), prec)
                if Vzz == 0.0:
                    eta = 0.0
                else:
                    eta = (Vxx - Vyy) / Vzz
                tmp = [Vxx, Vyy, Vzz, eta]
            else:
                Vzz = round(float(i.split()[3]), prec)
                Vyy = round(float(i.split()[2]), prec)
                Vxx = round(float(i.split()[1]), prec)
                if Vzz == 0.0:
                    eta = 0.0
                else:
                    eta = (Vyy - Vxx) / Vzz
                tmp = [Vxx, Vyy, Vzz, eta]
            efg_arr.append(tmp)
        efg_arr = np.array(efg_arr, dtype="float")
        return efg_arr

    @property
    def efg_raw_tensor(self):
        """Get raw electric field gradient tensor."""
        nions = self.nions
        for ii, i in enumerate(self.data):
            if "Electric field gradients (V/A^2)" in i:
                tmp = ii
        arr = self.data[tmp + 4 : tmp + 4 + nions]
        efg_arr = []
        for i in arr:
            line = i.split()
            tmp = [
                [line[1], line[4], line[5]],
                [line[4], line[2], line[6]],
                [line[5], line[6], line[3]],
            ]
            efg_arr.append(tmp)
        efg_arr = np.array(efg_arr, dtype="float")
        return efg_arr

    @property
    def quad_mom(self):
        """Get quadrupole momemnt."""
        nions = self.nions
        for ii, i in enumerate(self.data):
            if (
                "Q  : nuclear electric quadrupole moment in mb (millibarn)"
                in i
            ):
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
        """Get piezoelectric tensor."""
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
        Obtain elastic tensor and calculate related properties.

        For 3D and 2D cases.

        Args:
            outcar: OUTCAR file path

            vacuum: whether the structure has vaccum such as 2D materials
            for vacuum structures bulk and shear mod. Needs extra attention
            and elastic tensor are in Nm^-1 rather than GPa

        Returns:
              info: data for elastic tensor
              (in string and object representation),
              bulk, shear modulus, and phonon modes
        """
        ratio_c = 1.0
        if vacuum:
            ratio_c = 0.1 * float(
                abs(atoms.lattice_mat[2][2])
            )  # *(10**9)*(10**-10) #N/m unit
        KV = "na"
        GV = "na"
        # spin = "na"
        info = {}
        v = open(self.filename, "r")
        lines = v.read().splitlines()
        c = np.empty((6, 6), dtype=float)
        # TODO: Use regex to simplify
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
                    (c[0][0] + c[1][1] + c[2][2])
                    + 2 * (c[0][1] + c[1][2] + c[2][0])
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
        except Exception:
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
        """Initialize with filename."""
        with open(filename, "rb") as fp:

            def readData(dtype):
                """
                Read records from Fortran binary file.

                Convert to np.array of given dtype.
                """
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
        """Return the orbital derivative between states."""
        return self._cder_data

    @property
    def nbands(self):
        """Return the number of bands in the calculation."""
        return self._nbands

    @property
    def nkpoints(self):
        """Return the number of k-points in the calculation."""
        return self._nkpoints

    @property
    def nelect(self):
        """Return the number of electrons in the calculation."""
        return self._nelect

    def get_orbital_derivative_between_states(
        self, band_i, band_j, kpoint, spin, cart_dir
    ):
        """
        Return a valuebetween bands band_i and band_j.

        For k-point index, spin-channel and cartesian direction.

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
    Class for VASP Pseudowavefunction stored in WAVECAR.

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
        """Initialize the class."""
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
        except Exception:
            raise IOError("Failed to open %s" % self._fname)

        # read the basic information
        self.readWFHeader()
        # read the band information
        self.readWFBand()

        if self._lsoc:
            assert (
                self._nspin == 1
            ), "NSPIN = 1 for noncollinear version WAVECAR!"

    def isSocWfc(self):
        """Check if the WAVECAR is from an SOC calculation."""
        return True if self._lsoc else False

    def readWFHeader(self):
        """
        Read the system information from WAVECAR.

        It is written in the first two record.

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
        self._lattice_mat = dump[3:].reshape(
            (3, 3)
        )  # real space supercell basis
        self._Omega = np.linalg.det(
            self._lattice_mat
        )  # real space supercell volume
        self._Bcell = np.linalg.inv(
            self._lattice_mat
        ).T  # reciprocal space supercell volume

        # Minimum FFT grid size
        Anorm = np.linalg.norm(self._lattice_mat, axis=1)
        CUTOF = np.ceil(
            np.sqrt(self._encut / RYTOEV) / (TPI / (Anorm / AUTOA))
        )
        self._ngrid = np.array(2 * CUTOF + 1, dtype=int)

    def setWFPrec(self):
        """
        Set wavefunction coefficients precision.

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
        """Extract KS energies and Fermi occupations from WAVECAR."""
        self._nplws = np.zeros(self._nkpts, dtype=int)
        self._kvecs = np.zeros((self._nkpts, 3), dtype=float)
        self._energies = np.zeros(
            (self._nspin, self._nkpts, self._nbands), dtype=float
        )
        self._occs = np.zeros(
            (self._nspin, self._nkpts, self._nbands), dtype=float
        )

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
        Generate the G-vectors.

         satisfies the following relation
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
        """Read the planewave coefficients of specified KS states."""
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
        """Return the rec position for specified KS state."""
        self.checkIndex(ispin, ikpt, iband)

        rec = (
            2
            + (ispin - 1) * self._nkpts * (self._nbands + 1)
            + (ikpt - 1) * (self._nbands + 1)
            + iband
        )
        return rec

    def checkIndex(self, ispin, ikpt, iband):
        """Check if the index is valid."""
        assert 1 <= ispin <= self._nspin, "Invalid spin index!"
        assert 1 <= ikpt <= self._nkpts, "Invalid kpoint index!"
        assert 1 <= iband <= self._nbands, "Invalid band index!"


class Vasprun(object):
    """Construct vasprun.xml handling object."""

    def __init__(self, filename="vasprun.xml", data={}):
        """Intialize with filename, or optional parameters below."""
        self._filename = filename
        self._data = data
        self.ionic_steps = None
        self.electronic_steps = None
        self.input_parameters = None
        if self._data == {}:
            self.xml_to_dict()

    @classmethod
    def from_dict(self, d={}):
        """Construct class from a dictionary."""
        return Vasprun(filename=d["filename"], data=d["data"])

    def to_dict(self):
        """Convert class to a dictionary."""
        d = OrderedDict()
        d["filename"] = self._filename
        d["data"] = self._data
        return d

    def xml_to_dict(self):
        """Convert XML to dictionary."""
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
        """Get final energy."""
        return float(
            self.ionic_steps[-1]["scstep"][-1]["energy"]["i"][11]["#text"]
        )

    @property
    def nbands(self):
        """Get number of electronic bands."""
        return int(self.all_input_parameters["NBANDS"])

    @property
    def nkpoints(self):
        """Get number of kpoints."""
        return len(self.kpoints.kpts)

    @property
    def nspins(self):
        """Get total numnber of spins."""
        nspin = 1
        if self.is_spin_polarized:
            nspin = 2
        return nspin

    @property
    def efermi(self):
        """Get Fermi-energy."""
        return float(self.ionic_steps[-1]["dos"]["i"]["#text"])

    @property
    def num_atoms(self):
        """Get number of simulation atoms."""
        return int(self._data["modeling"]["atominfo"]["atoms"])

    @property
    def num_types(self):
        """Get number of atom types."""
        return int(self._data["modeling"]["atominfo"]["types"])

    @property
    def dielectric_loptics(self):
        """Get real and imag. dielectric function data."""
        tmp = self.ionic_steps[-1]["dielectricfunction"]["real"]["array"][
            "set"
        ]["r"]
        reals = []
        for i in range(len(tmp)):
            reals.append([float(j) for j in tmp[i].split()])

        tmp = self.ionic_steps[-1]["dielectricfunction"]["imag"]["array"][
            "set"
        ]["r"]
        imags = []
        for i in range(len(tmp)):
            imags.append([float(j) for j in tmp[i].split()])
        reals = np.array(reals)
        imags = np.array(imags)
        return reals, imags

    @property
    def avg_absorption_coefficient(self, max_axis=3):
        """Get average absoprtion coefficient. Used in solar-cell module."""
        eV_to_recip_cm = 1.0 / (
            physical_constants["Planck constant in eV s"][0]
            * speed_of_light
            * 1e2
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

    def phonon_data(self, fc_mass=True):
        """Get phonon data."""
        info = {}
        hessian = []
        # data = self._data
        for i in (self.ionic_steps[-1]["dynmat"]["varray"])[0]["v"]:
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
        if fc_mass:
            for i in range(natoms):
                for j in range(natoms):
                    force_constants[i, j] *= -np.sqrt(masses[i] * masses[j])
        phonon_eigenvals = np.array(
            self.ionic_steps[-1]["dynmat"]["v"]["#text"].split(),
            dtype="double",
        )
        eigvecs = np.array(
            [
                i.split()
                for i in (self.ionic_steps[-1]["dynmat"]["varray"][1]["v"])
            ],
            dtype="float",
        )
        phonon_eigenvectors = []
        for ev in eigvecs:
            phonon_eigenvectors.append(np.array(ev).reshape(natoms, 3))
        info["phonon_eigenvectors"] = phonon_eigenvectors
        info["phonon_eigenvalues"] = phonon_eigenvals
        info["masses"] = masses
        info["force_constants"] = force_constants
        return info

    @property
    def dfpt_data(self, fc_mass=True):
        """Get DFPT IBRION=8 LEPSILON lated data."""
        info = self.phonon_data(fc_mass=fc_mass)
        data = self._data
        natoms = self.num_atoms
        born_charges = []
        for n in range(natoms):
            born_charges.append(
                [
                    i.split()
                    for i in (data["modeling"]["calculation"]["array"]["set"])[
                        n
                    ]["v"]
                ]
            )
        born_charges = np.array(born_charges, dtype="double")
        epsilon = {}
        for i in data["modeling"]["calculation"]["varray"]:
            if "epsilon" in i["@name"]:
                epsilon[i["@name"]] = np.array(
                    [j.split() for j in i["v"]], dtype="float"
                )
        info["born_charges"] = born_charges
        info["epsilon"] = epsilon
        return info

    @property
    def get_dir_gap(self):
        """Get direct bandgap."""
        if not self.is_spin_polarized:
            spin_channels = 2
            # up_gap = "na"
            # dn_gap = "na"
            if self.is_spin_orbit:
                spin_channels = 1
            levels = int(
                float(self.all_input_parameters["NELECT"])
                / float(spin_channels)
            )
            ups = self.eigenvalues[0][:, :, 0][:, levels]
            dns = self.eigenvalues[0][:, :, 0][:, levels - 1]
            gap = min(ups - dns)
        else:
            tmp = np.concatenate(
                (self.eigenvalues[0][:, :, 0], self.eigenvalues[1][:, :, 0]),
                axis=1,
            )
            cat = np.sort(tmp, axis=1)
            nelect = int(float(self.all_input_parameters["NELECT"]))
            ups = cat[:, nelect]
            dns = cat[:, nelect - 1]
            gap = min(ups - dns)

        return gap

    @property
    def get_indir_gap(self):
        """Get indirect bandgap."""
        if not self.is_spin_polarized:
            spin_channels = 2
            # up_gap = "na"
            # dn_gap = "na"
            if self.is_spin_orbit:
                spin_channels = 1
            levels = int(
                float(self.all_input_parameters["NELECT"])
                / float(spin_channels)
            )
            gap = min(self.eigenvalues[0][:, :, 0][:, levels]) - max(
                self.eigenvalues[0][:, :, 0][:, levels - 1]
            )
            vbm = max(self.eigenvalues[0][:, :, 0][:, levels - 1])
            cbm = min(self.eigenvalues[0][:, :, 0][:, levels])

        if self.is_spin_polarized:
            tmp = np.concatenate(
                (self.eigenvalues[0][:, :, 0], self.eigenvalues[1][:, :, 0]),
                axis=1,
            )
            cat = np.sort(tmp, axis=1)
            nelect = int(float(self.all_input_parameters["NELECT"]))
            cbm = max(cat[:, nelect - 1])
            vbm = min(cat[:, nelect])
            gap = min(cat[:, nelect]) - max(cat[:, nelect - 1])
        return gap, cbm, vbm

    @property
    def fermi_velocities(self):
        """Get fermi velocities in m/s."""
        # TODO: check for other materials than graphene
        fermi_velocities = []
        fermi_k = []
        bands_cross_fermi = []
        h_bar = 6.582119569e-16  # reduced Planck const. eV s
        strt = self.all_structures[-1]
        lat = strt.lattice.reciprocal_lattice()
        kpoints_frac = self.kpoints.kpts
        kpoints_cart = [lat.cart_coords(i) for i in kpoints_frac]
        kpoints = np.array(kpoints_cart)
        for i, ii in enumerate(self.eigenvalues):
            for j, jj in enumerate(ii.T):
                for k, kk in enumerate(jj):
                    if max(kk) > self.efermi and min(kk) < self.efermi:
                        bands_cross_fermi.append(kk)
        for i in bands_cross_fermi:
            for j, jj in enumerate(i):
                if j < len(i) - 1:
                    if (i[j] < self.efermi < i[j + 1]) or (
                        i[j] > self.efermi > i[j + 1]
                    ):
                        # dk = np.sqrt(
                        #     (kpoints[j + 1][0] - kpoints[j][0]) ** 2
                        #     + (kpoints[j + 1][1] - kpoints[j][1]) ** 2
                        # )
                        dk = np.linalg.norm(kpoints[j + 1] - kpoints[j])
                        v_f = abs((i[j + 1] - i[j]) / (h_bar * dk))
                        fermi_velocities.append(v_f)
                        fermi_k.append(kpoints[j])
        # Convert to m/s
        # For graphene: ~0.85e6 m/s
        fermi_velocities = 1e-10 * np.array(fermi_velocities)
        return fermi_velocities, fermi_k, bands_cross_fermi

    @property
    def elements(self):
        """Get atom elements."""
        element_dat = self._data["modeling"]["atominfo"]["array"][0]["set"][
            "rc"
        ]
        if isinstance(element_dat, list):
            elements = [
                (element_dat[i]["c"][0]) for i in range(len(element_dat))
            ]
        elif isinstance(element_dat, dict):
            elements = [element_dat["c"][0]]
        else:
            raise ValueError("Unknown element type")
        if len(elements) != self.num_atoms:
            ValueError("Number of atoms is  not equal to number of elements")
        elements = [str(i) for i in elements]
        # print ('elements',elements)
        return elements

    def vrun_structure_to_atoms(self, s={}):
        """Convert structure to Atoms object."""
        tmp = s["crystal"]["varray"][0]["v"]
        lattice_mat = np.array([[float(j) for j in i.split()] for i in tmp])
        coord_info = s["varray"]["v"]
        # print ('coord_info',coord_info,type(coord_info))
        if isinstance(coord_info, list):
            frac_coords = np.array(
                [[float(j) for j in i.split()] for i in coord_info],
                dtype="float",
            )
        elif isinstance(coord_info, str):
            frac_coords = np.array(coord_info.split(), dtype="float")
            # print ('frac_coords',frac_coords)
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
        """Get all total energies."""
        energies = []
        for i in self.ionic_steps:
            en = float(i["energy"]["i"][1]["#text"])
            energies.append(en)
        return np.array(energies)

    @property
    def is_spin_polarized(self):
        """Check if the calculation is spin polarized."""
        if self.all_input_parameters["ISPIN"] == "2":
            return True
        else:
            return False

    @property
    def is_spin_orbit(self):
        """Check if the calculation is spin orbit."""
        if self.all_input_parameters["LSORBIT"] == "T":
            return True
        else:
            return False

    @property
    def all_structures(self):
        """Get all structures."""
        structs = []
        for i in self.ionic_steps:
            s = i["structure"]
            atoms = self.vrun_structure_to_atoms(s)
            structs.append(atoms)
        return structs

    @property
    def eigenvalues(self):
        """Get all eigenvalues."""
        nkpts = len(self.kpoints._kpoints)
        all_up_eigs = []
        all_dn_eigs = []
        if self.is_spin_polarized:
            for j in range(nkpts):
                eigs = np.array(
                    [
                        [float(jj) for jj in ii.split()]
                        for ii in (
                            self.ionic_steps[-1]["eigenvalues"]["array"][
                                "set"
                            ]["set"][0]
                        )["set"][j]["r"]
                    ]
                )
                all_up_eigs.append(eigs)
            for j in range(nkpts):
                eigs = np.array(
                    [
                        [float(jj) for jj in ii.split()]
                        for ii in (
                            self.ionic_steps[-1]["eigenvalues"]["array"][
                                "set"
                            ]["set"][1]
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
                            self.ionic_steps[-1]["eigenvalues"]["array"][
                                "set"
                            ]["set"]
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
        """Get all forces."""
        forces = []
        for m in self.ionic_steps:
            force = np.array(
                [[float(j) for j in i.split()] for i in m["varray"][0]["v"]]
            )

            forces.append(force)
        return np.array(forces)

    @property
    def all_stresses(self):
        """Get all stresses."""
        stresses = []
        for m in self.ionic_steps:
            stress = np.array(
                [[float(j) for j in i.split()] for i in m["varray"][1]["v"]]
            )

            stresses.append(stress)
        return np.array(stresses)

    @property
    def all_input_parameters(self):
        """Get all explicit input parameters. Need to add a few more."""
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
        """Get Kpoints."""
        kplist = np.array(
            [
                [float(j) for j in i.split()]
                for i in self._data["modeling"]["kpoints"]["varray"][0]["v"]
            ]
        )
        kpwt = np.array(
            [
                float(i)
                for i in self._data["modeling"]["kpoints"]["varray"][1]["v"]
            ]
        )
        return Kpoints(kpoints=kplist, kpoints_weights=kpwt)

    def get_bandstructure(
        self,
        E_low=-4,
        E_high=4,
        spin=0,
        zero_efermi=True,
        kpoints_file_path="KPOINTS",
        plot=False,
    ):
        """Get electronic bandstructure plot."""
        try:
            kp_labels = []
            kp_labels_points = []
            f = open(kpoints_file_path, "r")
            lines = f.read().splitlines()
            f.close()
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

        except Exception:
            print("No K-points file found, still proceeding")
            pass

        tmp = 0.0
        info = {}
        info["efermi"] = float(self.efermi)
        if zero_efermi:
            tmp = float(self.efermi)

        spin_up_bands_x = []
        spin_up_bands_y = []
        spin_down_bands_x = []
        spin_down_bands_y = []
        for i, ii in enumerate(self.eigenvalues[spin][:, :, 0].T - tmp):
            # plt.plot(ii, color="r")
            spin_up_bands_x.append([np.arange(0, len(ii))])
            spin_up_bands_y.append([ii])
        if self.is_spin_polarized:
            for i, ii in enumerate(self.eigenvalues[1][:, :, 0].T - tmp):
                # plt.plot(ii, color="b")
                spin_down_bands_x.append([np.arange(0, len(ii))])
                spin_down_bands_y.append([ii])

        info["spin_up_bands_x"] = spin_up_bands_x
        info["spin_up_bands_y"] = spin_up_bands_y
        info["spin_down_bands_x"] = spin_down_bands_x
        info["spin_down_bands_y"] = spin_down_bands_y

        info["kp_labels_points"] = list(kp_labels_points)
        info["kp_labels"] = list(kp_labels)
        if plot:
            for i, j in zip(info["spin_up_bands_x"], info["spin_up_bands_y"]):
                plt.plot(
                    np.array(i).flatten(), np.array(j).flatten(), color="b"
                )

            if self.is_spin_polarized:
                for i, j in zip(
                    info["spin_down_bands_x"], info["spin_down_bands_y"]
                ):
                    plt.plot(
                        np.array(i).flatten(), np.array(j).flatten(), color="r"
                    )

            plt.ylim([E_low, E_high])
            plt.xticks(kp_labels_points, kp_labels)
            plt.xlim([0, len(self.kpoints._kpoints)])
            plt.xlabel(r"$\mathrm{Wave\ Vector}$")
            ylabel = (
                r"$\mathrm{E\ -\ E_f\ (eV)}$"
                if zero_efermi
                else r"$\mathrm{Energy\ (eV)}$"
            )
            plt.ylabel(ylabel)
        return info

    @property
    def total_dos(self):
        """Get total density of states."""
        energies = []
        spin_up = []
        spin_dn = []
        spin_up_data = np.array(
            [
                [float(j) for j in i.split()]
                for i in self.ionic_steps[-1]["dos"]["total"]["array"]["set"][
                    "set"
                ][0]["r"]
            ]
        )
        energies = spin_up_data[:, 0]
        spin_up = spin_up_data[:, 1]
        if self.is_spin_polarized:
            spin_dn = []
            spin_dn_data = np.array(
                [
                    [float(j) for j in i.split()]
                    for i in self.ionic_steps[-1]["dos"]["total"]["array"][
                        "set"
                    ]["set"][1]["r"]
                ]
            )
            spin_dn = -1 * spin_dn_data[:, 1]
        return energies, spin_up, spin_dn

    @property
    def partial_dos_spdf(self):
        """Get partial density of states."""
        info = rec_dict()
        natoms = self.num_atoms
        nspin = self.nspins
        pdos_keys = self.ionic_steps[-1]["dos"]["partial"]["array"]["field"]
        steps_dat = self.ionic_steps[-1]["dos"]["partial"]["array"]["set"][
            "set"
        ]
        if isinstance(steps_dat, list):
            for atom in range(natoms):
                for spin in range(nspin):
                    for k, key in enumerate(pdos_keys):
                        vals = np.array(
                            [
                                ii.split()
                                for ii in (
                                    self.ionic_steps[-1]["dos"]["partial"][
                                        "array"
                                    ]["set"]["set"][atom]["set"][spin]["r"]
                                )
                            ],
                            dtype="float",
                        )

                        info[spin][atom][key] = vals[:, k]
        elif isinstance(steps_dat, dict):
            atom = 0
            for spin in range(nspin):
                for k, key in enumerate(pdos_keys):
                    vals = np.array(
                        [
                            ii.split()
                            for ii in (
                                self.ionic_steps[-1]["dos"]["partial"][
                                    "array"
                                ]["set"]["set"]["set"][spin]["r"]
                            )
                        ],
                        dtype="float",
                    )

                    info[spin][atom][key] = vals[:, k]
        else:
            raise ValueError("Bug in PDOS parser.")

        return info

    @property
    def projected_spins_kpoints_bands(self):
        """Use for spin, kpoint and band projected bandstructure plots."""
        info = rec_dict()
        nspin = self.nspins
        nkpoints = self.nkpoints
        nbands = self.nbands
        for spin in range(nspin):
            for kpoint in range(nkpoints):
                for nb in range(nbands):
                    vals = [
                        float(ii)
                        for ii in (
                            self.ionic_steps[-1]["projected"]["eigenvalues"][
                                "array"
                            ]["set"]["set"][spin]["set"][kpoint]["r"][nb]
                        ).split()
                    ]
                    info[spin][kpoint][nb] = vals
        return info

    @property
    def projected_atoms_spins_kpoints_bands(self):
        """Use for atom,spin,kpoint and band projected bandstructures."""
        info = rec_dict()
        # orbitals = self.ionic_steps[-1]["projected"]["array"]["field"]
        dimensions = [
            list(i.values())[1]
            for i in self.ionic_steps[-1]["projected"]["array"]["dimension"]
        ]
        natoms = self.num_atoms
        nkpoints = self.nkpoints
        nbands = self.nbands
        nspin = self.nspins
        for atom in range(natoms):
            for spin in range(nspin):
                for kpoint in range(nkpoints):
                    for band in range(nbands):
                        for orbital in dimensions:
                            val = (
                                (
                                    self.ionic_steps[-1]["projected"]["array"][
                                        "set"
                                    ]["set"][spin]
                                )["set"][kpoint]["set"][band]
                            )["r"][atom].split()
                            val = [float(v) for v in val]
                            info[atom][spin][kpoint][band][orbital] = np.array(
                                val
                            )
        return info

    def get_spdf_dos(self, plot=False):
        """Get spdf resolved partial density of states."""
        info = {}
        s = ["s"]
        p = ["px", "py", "pz"]
        d = ["dxy", "dyz", "dxz", "dz2", "dx2"]
        f = ["f-3", "f-2", "f-1", "f0", "f1", "f2", "f3"]
        # spin = 0
        spin_pol = self.is_spin_polarized
        num_atoms = self.all_structures[-1].num_atoms
        pdos = self.partial_dos_spdf  # spin,atom,spdf
        energy = pdos[0][0]["energy"] - self.efermi
        spin_up_s = np.zeros(len(energy))
        spin_up_p = np.zeros(len(energy))
        spin_up_d = np.zeros(len(energy))
        has_f_elements = False
        if "f1" in pdos[0][0].keys():
            has_f_elements = True
            spin_up_f = np.zeros(len(energy))
        for i in range(num_atoms):
            for j in s:
                if isinstance(pdos[0][i][j], np.ndarray):
                    spin_up_s += pdos[0][i][j]
            for j in p:
                if isinstance(pdos[0][i][j], np.ndarray):
                    spin_up_p += pdos[0][i][j]
            for j in d:
                if isinstance(pdos[0][i][j], np.ndarray):
                    spin_up_d += pdos[0][i][j]
            if has_f_elements:
                for j in f:
                    if isinstance(pdos[0][i][j], np.ndarray):
                        spin_up_f += pdos[0][i][j]
        info["spin_up_s"] = spin_up_s
        info["spin_up_p"] = spin_up_p
        info["spin_up_d"] = spin_up_d
        info["energy"] = energy
        if has_f_elements:
            info["spin_up_f"] = spin_up_f
        if spin_pol:
            # spin = 1
            spin_down_s = np.zeros(len(energy))
            spin_down_p = np.zeros(len(energy))
            spin_down_d = np.zeros(len(energy))
            if has_f_elements:
                spin_down_f = np.zeros(len(energy))
            for i in range(num_atoms):
                for j in s:
                    if isinstance(pdos[0][i][j], np.ndarray):
                        spin_down_s += pdos[1][i][j]
                for j in p:
                    if isinstance(pdos[0][i][j], np.ndarray):
                        spin_down_p += pdos[1][i][j]
                for j in d:
                    if isinstance(pdos[0][i][j], np.ndarray):
                        spin_down_d += pdos[1][i][j]
                if has_f_elements:
                    for j in f:
                        if isinstance(pdos[0][i][j], np.ndarray):
                            spin_down_f += pdos[1][i][j]

            info["spin_down_s"] = -1 * spin_down_s
            info["spin_down_p"] = -1 * spin_down_p
            info["spin_down_d"] = -1 * spin_down_d
            info["has_f_elements"] = str(has_f_elements)
            if has_f_elements:
                info["spin_down_f"] = -1 * spin_down_f
            if plot:
                plt.plot(
                    info["energy"], info["spin_up_s"], color="red", label="s"
                )

                plt.plot(
                    info["energy"], info["spin_up_p"], color="green", label="p"
                )

                plt.plot(
                    info["energy"], info["spin_up_d"], color="blue", label="d"
                )
                if has_f_elements:
                    plt.plot(
                        info["energy"],
                        info["spin_up_f"],
                        color="black",
                        label="f",
                    )
                if spin_pol:
                    plt.plot(info["energy"], info["spin_down_s"], color="red")
                    plt.plot(
                        info["energy"], info["spin_down_p"], color="green"
                    )
                    plt.plot(info["energy"], info["spin_down_d"], color="blue")
                if has_f_elements:
                    plt.plot(
                        info["energy"], info["spin_down_f"], color="black"
                    )
                plt.xlim([-5, 10])
                plt.legend()
        return info

    def get_atom_resolved_dos(self, plot=False):
        """Get atom resolved density of states."""
        # spin_pol = self.is_spin_polarized
        atoms = self.all_structures[-1]
        # num_atoms = atoms.num_atoms
        elements = atoms.elements
        unique_elements = atoms.uniq_species
        pdos = self.partial_dos_spdf  # spin,atom,spdf
        energy = pdos[0][0]["energy"] - self.efermi
        element_dict = recast_array_on_uniq_array_elements(
            unique_elements, elements
        )
        valid_keys = []
        info = {}
        info["spin_up_info"] = {}
        info["spin_down_info"] = {}
        info["energy"] = energy
        for i in pdos[0][0].keys():
            if "energ" not in i:
                valid_keys.append(i)
        # print (valid_keys)
        spin_up_info = {}
        for i, j in element_dict.items():
            spin_up_info[i] = np.zeros(len(energy))

        for i, j in element_dict.items():
            for atom in j:
                for k in valid_keys:
                    spin_up_info[i] += pdos[0][atom][k]
        info["spin_up_info"] = spin_up_info
        if self.is_spin_polarized:
            spin_down_info = {}
            for i, j in element_dict.items():
                spin_down_info[i] = np.zeros(len(energy))

            for i, j in element_dict.items():
                for atom in j:
                    for k in valid_keys:
                        spin_down_info[i] += -1 * pdos[0][atom][k]
            info["spin_down_info"] = spin_down_info
            if plot:
                for i, j in info.items():
                    if "spin" in i:
                        for m, n in j.items():

                            if "up" in i:
                                plt.plot(info["energy"], n, label=m)
                            if "down" in i:
                                plt.plot(info["energy"], n)
                plt.legend()
        return info


def parse_raman_dat(
    vasp_raman_path="RAMANDIR-bulk@JVASP-1002_mp-149/vasp_raman.dat",
):
    """
    Parse vasp_raman.dat .

    generated by https://github.com/raman-sc/VASP
    """
    f = open(vasp_raman_path, "r")
    lines = f.read().splitlines()
    f.close()
    info = {}
    freqs = []
    activity = []
    alpha = []
    beta2 = []
    for i in lines:
        if "#" not in i:
            tmp = i.split()
            freqs.append(float(tmp[1]))
            alpha.append(float(tmp[2]))
            beta2.append(float(tmp[3]))
            activity.append(float(tmp[4]))
    freqs = np.array(freqs)
    activity = np.array(activity)
    indices = np.arange(0, len(activity) - 1)
    try:
        indices = ss.find_peaks_cwt(activity, np.arange(1, 5))
    except Exception:
        print("Cannot use peak finding module", vasp_raman_path)
        pass
    info["freqs"] = freqs
    info["activity"] = activity
    info["alpha"] = alpha
    info["beta2"] = beta2
    info["indices"] = indices
    return info


"""
kp='/users/knc6/Software/Devs/jarvis/jarvis/examples/vasp/SiOptb88/MAIN-RELAX-bulk@mp_149/KPOINT'
kpt=Kpoints(filename=kp)
print (kpt)
"""
