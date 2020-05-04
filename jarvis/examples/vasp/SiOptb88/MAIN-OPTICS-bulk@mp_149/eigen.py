from pymatgen.io.vasp.outputs import Vasprun
v=Vasprun('vasprun.xml', occu_tol=1e-2)
print v.eigenvalue_band_properties[0]
