import os
from jarvis.io.phonopy.outputs import bandstructure_plot,total_dos
totdos = os.path.join(os.path.dirname(__file__),'total_dos.dat')
band = os.path.join(os.path.dirname(__file__),'band.yaml')
def test_outputs():
  #JVASP-1002
  x,y=total_dos(tot_dos=totdos)
  a,b,c,d=bandstructure_plot(band_yaml=band)
 
