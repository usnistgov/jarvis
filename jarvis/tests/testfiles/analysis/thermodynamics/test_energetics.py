from jarvis.analysis.thermodynamics.energetics import form_enp,get_twod_defect_energy
import os
from jarvis.io.vasp.outputs import Vasprun
tmp_xml=os.path.join(os.path.dirname(__file__), "JVASP-667_C_C_c.xml")
vrun=Vasprun(tmp_xml)
def test_get_twod_defect_energy():
  Ef=get_twod_defect_energy(vrun=vrun,jid="JVASP-667",atom="C")
  print (Ef)
def test_form_enp():
   atoms=vrun.all_structures[-1]
   total_energy=vrun.final_energy
   Ef=form_enp(atoms=atoms,total_energy=total_energy)
   print (Ef)
#test_get_twod_defect_energy()
#test_form_enp()

