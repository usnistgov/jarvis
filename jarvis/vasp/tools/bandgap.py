from pymatgen.io.vasp.outputs import Vasprun
v=Vasprun('vasprun.xml',occu_tol=0.001)
print ("bandgap",v.eigenvalue_band_properties[0])
