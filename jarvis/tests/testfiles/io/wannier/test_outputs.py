from jarvis.io.wannier.outputs import (
    WannierHam,
    Wannier90wout,
    Wannier90eig,
    get_projectors_for_formula,
    get_orbitals,
)
import os
import tempfile
from jarvis.core.kpoints import generate_kgrid
from jarvis.io.vasp.inputs import Poscar
import matplotlib.pyplot as plt

plt.switch_backend("agg")

new_file, filename = tempfile.mkstemp()
atoms = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "analysis",
        "structure",
        "POSCAR",
    )
).atoms
wann_soc_win_hr = os.path.join(os.path.dirname(__file__), "wannier90_hr.dat")
wann_wout = os.path.join(os.path.dirname(__file__), "wannier90.wout")
wann_eig = os.path.join(os.path.dirname(__file__), "wannier90.eig")
soc_scfband_vrun = os.path.join(
    os.path.dirname(__file__), "vasprun.xml"
)  # for JVASP-1067


def test_outputs_bi2se3():
    pp = get_projectors_for_formula(formula_dict={"Cr": 1, "I": 3})[0]
    orb = get_orbitals()[0]
    x = get_orbitals(
        projection_info=[["Cr", 4, ["s", "d"]], ["I", 3, ["s", "p"]]],
        desired_orbitals=[["Cr", "d"]],
    )
    print(x)
    x = get_orbitals(
        projection_info=[["Cr", 4, ["s", "d"]], ["I", 3, ["s", "p"]]],
        desired_orbitals=[["Cr", "s"]],
    )
    print(x)
    w = WannierHam(filename=wann_soc_win_hr)
    new_file, filename = tempfile.mkstemp()
    comp = w.compare_dft_wann(
        vasprun_path=soc_scfband_vrun, plot=True, filename=filename + ".png"
    )
    maxdiff = comp["maxdiff"]
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
    new_file, filename = tempfile.mkstemp()
    # print(pp, orb)
    w.get_bandstructure_plot(atoms=atoms, filename=filename)
    # print (x,pp)
    w.find_nodes(nk1=1,nk2=1,nk3=1)
    w.fermi_surf_2d(nk1=1,nk2=1)
    w.chern_number_simple()
    # print (round(dos[75],3))
    assert (
        round(maxdiff, 2),
        round(dos[75], 3),
        pp[0][0],
        orb,
        w.nwan,
        big.nwan,
    ) == (0.12, 3.02, "Bi", 1, 40, 320,)


def test_wann_cent():
    centers = Wannier90wout(wout_path=wann_wout).give_wannier_centers()
    # print (centers, len(centers))
    assert (len(centers)) == (40)

def wann_eig():
    eigs = Wannier90eig(wann_eig)
    eigs.give_wannier_eigs()
    eigs.neigs()
    eigs.nk()

# test_outputs_cri3()
# test_wann_cent()
