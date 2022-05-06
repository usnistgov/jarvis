"""Module for processing elastic tensor."""
# Reference: https://doi.org/10.1103/PhysRevB.98.014107
import numpy as np
from collections import OrderedDict


class ElasticTensor(object):
    """Module for processing elastic tensor."""

    def __init__(self, et_tensor=[]):
        """Initialize class."""
        self.et_tensor = et_tensor

    @property
    def voigt_modulus(self):
        """Get Voigt modulus."""
        c = self.et_tensor
        Kv = float(
            (c[0][0] + c[1][1] + c[2][2]) + 2 * (c[0][1] + c[1][2] + c[2][0])
        ) / float(9)
        Gv = float(
            (c[0][0] + c[1][1] + c[2][2])
            - (c[0][1] + c[1][2] + c[2][0])
            + 3 * (c[3][3] + c[4][4] + c[5][5])
        ) / float(15)
        return [Kv, Gv]

    @property
    def compliance_tensor(self):
        """Get compliance."""
        return np.linalg.inv(self.et_tensor)

    @property
    def reuss_modulus(self):
        """Get Reuss modulus."""
        c = self.compliance_tensor
        Kr = 1 / float(
            (c[0][0] + c[1][1] + c[2][2]) + 2 * (c[0][1] + c[1][2] + c[2][0])
        )
        Gr = 15 / (
            4 * (c[0][0] + c[1][1] + c[2][2])
            - 4 * (c[0][1] + c[1][2] + c[2][0])
            + 3 * (c[3][3] + c[4][4] + c[5][5])
        )
        return [Kr, Gr]

    def long_velocity(self, atoms=None):
        """Longitudinal velocity using Navier equation."""
        # density = atoms.density
        weight = float(atoms.composition.weight)
        volume = atoms.volume
        mass_density = 1.6605e3 * weight / volume
        avg_mod = self.average_modulus
        k_vrh = avg_mod[0]
        g_vrh = avg_mod[1]
        # 1e9:= GPa to Pascal (kg/ms^2)
        vel = np.sqrt(1e9 * (k_vrh + 4.0 / 3.0 * g_vrh) / mass_density)
        return vel

    def trans_velocity(self, atoms=None):
        """Transverse velocity."""
        avg_mod = self.average_modulus
        g_vrh = avg_mod[1]
        volume = atoms.volume
        weight = float(atoms.composition.weight)
        mass_density = 1.6605e3 * weight / volume
        vel = np.sqrt(1e9 * g_vrh / mass_density)
        return vel

    def velocity_average(self, atoms=None):
        """Average velocity."""
        vt = self.trans_velocity(atoms=atoms)
        vl = self.long_velocity(atoms=atoms)
        return 1.0 / (
            np.cbrt(
                (1.0 / 3.0) * (2.0 / (vt * vt * vt) + 1.0 / (vl * vl * vl))
            )
        )

    def debye_temperature(self, atoms=None):
        """Debye temperature."""
        const = 1.05457e-34 / 1.38065e-23  # (h/kb)
        v0 = atoms.volume * 1e-30 / atoms.num_atoms
        vl = self.long_velocity(atoms=atoms)
        vt = self.trans_velocity(atoms=atoms)
        vm = 3 ** (1.0 / 3.0) * (1 / vl ** 3 + 2 / vt ** 3) ** (-1.0 / 3.0)
        theta = const * vm * (6 * np.pi ** 2 / v0) ** (1.0 / 3.0)
        return theta

    @property
    def average_modulus(self):
        """Get average modulus."""
        return (
            np.array(self.voigt_modulus) + np.array(self.reuss_modulus)
        ) / 2

    @property
    def pugh_ratio_voigt(self):
        """Get Voigt Pugh ratio."""
        Kv, Gv = self.voigt_modulus
        return Gv / Kv

    @property
    def pettifor_criteria(self):
        """Get Pettifor criteria."""
        c = self.et_tensor
        return c[0][1] - c[3][3]

    @property
    def is_brittle(self):
        """Check if brittle material."""
        return self.pugh_ratio_voigt > 0.571 and self.pettifor_criteria < 0

    @property
    def is_ductile(self):
        """Check if ductile material."""
        return self.pugh_ratio_voigt < 0.571 and self.pettifor_criteria > 0

    @property
    def melting_temperature_metals(self):
        """Get crude Melting temp. estimate."""
        # https://doi.org/10.1016/0036-9748(84)90267-9
        avg_mod = self.average_modulus
        k_vrh = avg_mod[0]
        return 607 + 9.3 * k_vrh

    @property
    def cauchy_pressure(self):
        """Get Cauchy pressure."""
        # >0 ionic bonding
        # <0 covalent bonding
        c = self.et_tensor
        return c[0][1] - c[3][3]

    @property
    def poisson_ratio(self):
        """Get poisson's ratio."""
        k, g = self.average_modulus
        return (3 * k - 2 * g) / (6 * k + 2 * g)

    @property
    def universal_ansiotropy_ratio(self):
        """Get universal ansiotropy ratio."""
        Kv, Gv = self.voigt_modulus
        Kr, Gr = self.reuss_modulus
        return 5 * (Gv / Gr) + (Kv / Kr) - 6

    @property
    def youngs_modulus(self):
        """Get Youngs modulus."""
        k, g = self.average_modulus
        return 9e9 * k * g / (3 * k + g)

    def to_dict(self):
        """Get dictionary representation."""
        d = OrderedDict()
        d["voigt_bulk_modulus"] = self.voigt_modulus[0]
        d["voigt_shear_modulus"] = self.voigt_modulus[1]
        d["reuss_bulk_modulus"] = self.reuss_modulus[0]
        d["reuss_shear_modulus"] = self.reuss_modulus[1]
        d["poisson_ratio"] = self.poisson_ratio
        d["youngs_modulus"] = self.youngs_modulus
        d["universal_ansiotropy_ratio"] = self.universal_ansiotropy_ratio
        if not isinstance(self.et_tensor, list):
            et_tensor = self.et_tensor.tolist()
        else:
            et_tensor = self.et_tensor
        d["raw_et_tensor"] = et_tensor
        return d


"""
from jarvis.io.vasp.outputs import Vasprun,Outcar
o=Outcar('../../examples/vasp/SiOptb88/SiOptb88/MAIN-ELASTIC-bulk@mp_149/OUTCAR')
print (o.elastic_props())
p=Atoms.from_poscar('../../examples/vasp/SiOptb88/SiOptb88/MAIN-ELASTIC-bulk@mp_149/POSCAR')
et=ElasticTensor(o.elastic_props()['cij'])
#print (et.voigt_modulus) #[87.26666666666667, 63.28]
#print (et.reuss_modulus) #[87.26666666666665, 60.24546397096941]
#print (et.average_modulus) #[87.26666667 61.76273199]
#print (et.poisson_ratio) #0.21367500388646996
#print (et.universal_ansiotropy_ratio) #0.21367500388646996
#print ('j_vel',(1e9*(et.average_modulus[0])/p.density)**.5)
#print ('j_vel',(1e9*
(et.average_modulus[0]+4/3*et.average_modulus[1])/p.density)**.5)
k,g=et.average_modulus
mass_density=p.density
#print ((1e9 * (k + 4. / 3. * g) / mass_density/1.6605e3 ) ** 0.5)
print (et.to_dict())
from pymatgen.analysis.elasticity.elastic import ElasticTensor
et=ElasticTensor.from_voigt(o.elastic_props()['cij'])
from pymatgen.core.structure import Structure
pmg=Structure.from_file('../../examples/vasp/
SiOptb88/SiOptb88/MAIN-ELASTIC-bulk@mp_149/POSCAR')
#print (et.k_vrh,et.g_vrh)
#print (et.long_v(pmg))
#print ('density_j,density_p',p.density,pmg.density)
"""
