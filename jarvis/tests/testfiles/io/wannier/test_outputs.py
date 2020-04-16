from jarvis.io.wannier.outputs import (
    WannierHam,
    Wannier90wout,
    get_projectors_for_formula,
    get_orbitals,
)
import os
from jarvis.core.kpoints import generate_kgrid

wann_soc_win_hr = os.path.join(os.path.dirname(__file__), "wannier90_hr.dat")
wann_wout = os.path.join(os.path.dirname(__file__), "wannier90.wout")
soc_scfband_vrun = os.path.join(
    os.path.dirname(__file__), "vasprun.xml"
)  # for JVASP-1067


def test_outputs():
    w = WannierHam(filename=wann_soc_win_hr)
    maxdiff = w.compare_dft_wann(vasprun_path=soc_scfband_vrun, plot=False)["maxdiff"]
    info = w.to_dict()
    dd = WannierHam.from_dict(info)

    pp = get_projectors_for_formula()[0][0]
    # print("getorbs", pp)
    kpoints = generate_kgrid([5, 5, 5])
    energies, dos, pdos = w.dos(kpoints)
    orb = get_orbitals()[0]
    big = w.generate_supercell([2, 2, 2])

    pp = get_projectors_for_formula()
    x = get_orbitals()
    print(pp, orb)

    # print (x,pp)

    # print (round(dos[75],3))
    assert (round(maxdiff, 2), round(dos[75], 3), pp[0][0], orb, w.nwan, big.nwan) == (
        0.12,
        2.881,
        "Bi",
        1,
        40,
        320,
    )


def test_wann_cent():
    centers = Wannier90wout(wout_path=wann_wout).give_wannier_centers()
    # print (centers, len(centers))
    assert (len(centers)) == (40)

test_outputs()
# test_wann_cent()
