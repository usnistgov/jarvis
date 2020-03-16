from jarvis.analysis.spillage.spillage import  Spillage
import os

wf_noso = os.path.join('.','WAVECAR.nosoc')
wf_so = os.path.join('.','WAVECAR.soc')


def test_spillage():
    spl = Spillage(wf_noso=wf_noso,wf_so=wf_so)
    info = spl.overlap_so_spinpol()
    spillage = round(info['spillage'],2)
    #print (spillage)
    assert (spillage)==(1.36)
#test_spillage()

