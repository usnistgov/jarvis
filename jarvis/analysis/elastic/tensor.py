"""Module for processing elastic tensor."""
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

    @property
    def average_modulus(self):
        """Get average modulus."""
        return (
            np.array(self.voigt_modulus) + np.array(self.reuss_modulus)
        ) / 2

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
        d["raw_et_tensor"] = self.et_tensor
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
