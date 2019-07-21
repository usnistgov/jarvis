from jarvis.spillage.vaspwfc import vaspwfc
import os

wf_so = os.path.join(os.path.dirname(__file__), "..", "spillage", "WAVECAR.soc")

def test_vaspwfc():
    so = vaspwfc(fnm=wf_so, lsorbit=True)
    so.printWF()
    so_k, so_bands = so.readWFBand()
    print (so_bands)
