from jarvis.io.wannier.outputs import WannierHam, Wannier90wout
import os

wann_soc_win_hr = os.path.join(os.path.dirname(__file__), "wannier90_hr.dat")
wann_wout = os.path.join(os.path.dirname(__file__), "wannier90.wout")
soc_scfband_vrun = os.path.join(
    os.path.dirname(__file__), "vasprun.xml"
)  # for JVASP-1067


def test_outputs():
    w = WannierHam(filename=wann_soc_win_hr)
    maxdiff = w.compare_dft_wann(vasprun_path=soc_scfband_vrun, plot=False)
    assert (round(maxdiff, 2)) == (0.12)


def test_wann_cent():
    centers = Wannier90wout(wout_path=wann_wout).give_wannier_centers()
    # print (centers, len(centers))
    assert (len(centers)) == (40)


# test_wann_cent()
