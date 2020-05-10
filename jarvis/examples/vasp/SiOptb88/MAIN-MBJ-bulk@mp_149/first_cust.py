"""
from pymatgen.io.vasp.inputs import Incar, Poscar, VaspInput,Potcar, Kpoints
import os,shutil
from custodian.vasp.jobs import VaspJob
from custodian.vasp.handlers import VaspErrorHandler, UnconvergedErrorHandler,MeshSymmetryErrorHandler, NonConvergingErrorHandler, PotimErrorHandler
from custodian.vasp.validators import VasprunXMLValidator
from custodian.custodian import Custodian
inc=Incar.from_file("INCAR")
pot=Potcar.from_file("POTCAR")
pos=Poscar.from_file("POSCAR")
kp=Kpoints.from_file("KPOINTS")
shutil.copy2('/users/knc6/bin/vdw_kernel.bindat','./')
vinput = VaspInput.from_directory(".")
job=VaspJob(['mpirun', '-np', '16', '/users/knc6/VASP/vasp54/src/vasp.5.4.1/bin/vasp_std'], final=False, backup=False)
handlers = [VaspErrorHandler(), MeshSymmetryErrorHandler(),UnconvergedErrorHandler(), NonConvergingErrorHandler(),PotimErrorHandler()]
validators = [VasprunXMLValidator()]
c = Custodian(handlers, [job],max_errors=5,validators=validators)
c.run()
"""
