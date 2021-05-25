from jarvis.analysis.phonon.force_constants import qpoint,read_fc
import os
fc_file = os.path.join(os.path.dirname(__file__), "FORCE_CONSTANTS")

def test_fc():
   fc=read_fc(fc_file)
   assert (fc[0][0][0][0])==12.960974735
   qp=qpoint(force_constant=fc,qpt=[0, 0, 0])
   assert (qp[0][0][0][0])==12.960974735
#test_fc()   
