from jarvis.io.vasp.outputs import (
    Vasprun,
    Oszicar,
    Wavecar,
    Waveder,
    Chgcar,
    Locpot,
    Outcar,
    parse_raman_dat,
)
import numpy as np
import os
from jarvis.analysis.phonon.ir import ir_intensity
import matplotlib.pyplot as plt

plt.switch_backend("agg")


import tarfile

example_fold_tgz = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88.tgz",
)


example_fold = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
)

if not os.path.isdir(example_fold):
    tar = tarfile.open(example_fold_tgz)
    tar.extractall(example_fold)
    tar.close()


vrun = Vasprun(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "SiOptb88",
        "MAIN-RELAX-bulk@mp_149",
        "vasprun.xml",
    )
)
band_vrun = Vasprun(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "SiOptb88",
        "MAIN-BAND-bulk@mp_149",
        "vasprun.xml",
    )
)
non_spinpol_vrun = Vasprun(
    filename=os.path.join(
        os.path.dirname(__file__),
        "vasprun.xml.JVASP-23436",
    )
)
vasp544_mbj_optics_vrun = Vasprun(
    filename=os.path.join(
        os.path.dirname(__file__),
        "vasprun.xml.JVASP-97577",
    )
)
opt_vrun = Vasprun(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "SiOptb88",
        "MAIN-OPTICS-bulk@mp_149",
        "vasprun.xml",
    )
)
band_kp = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
    "SiOptb88",
    "MAIN-BAND-bulk@mp_149",
    "KPOINTS",
)


loc = Locpot(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "SiOptb88",
        "MAIN-RELAX-bulk@mp_149",
        "LOCPOT",
    )
)


chg = Chgcar(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "SiOptb88",
        "MAIN-RELAX-bulk@mp_149",
        "CHGCAR",
    )
)
osz = Oszicar(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "SiOptb88",
        "MAIN-RELAX-bulk@mp_149",
        "OSZICAR",
    )
)
out = Outcar(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "SiOptb88",
        "MAIN-ELASTIC-bulk@mp_149",
        "OUTCAR",
    )
)
opt_out = Outcar(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "SiOptb88",
        "MAIN-OPTICS-bulk@mp_149",
        "OUTCAR",
    )
)
# TODO:
# wder = Waveder(
#    os.path.join(
#        os.path.dirname(__file__),
#        "..",
#        "..",
#        "..",
#        "..",
#        "examples",
#        "vasp",
#        "SiOptb88",
#        "SiOptb88",
#        "MAIN-OPTICS-bulk@mp_149",
#        "WAVEDER",
#    )
# )
wf_noso = Wavecar(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "analysis",
        "topological",
        "WAVECAR.nosoc",
    ),
    lsorbit=False,
)
wf_so = Wavecar(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "analysis",
        "topological",
        "WAVECAR.soc",
    ),
    lsorbit=True,
)


def test_chgcar():
    # print (chg.is_spin_polarized(), chg.is_spin_orbit(), np.array(chg.chg).shape)
    assert (
        chg.is_spin_polarized(),
        chg.is_spin_orbit(),
        np.array(chg.chg).shape,
    ) == (
        True,
        False,
        (2, 56, 56, 56),
    )
    td = chg.to_dict()
    fd = Chgcar.from_dict(td)
    x = chg.modify_grid()
    cmd = "rm New_CHGCAR"
    os.system(cmd)


def test_locpot():
    # print (chg.is_spin_polarized(), chg.is_spin_orbit(), np.array(chg.chg).shape)
    assert (
        loc.is_spin_polarized(),
        loc.is_spin_orbit(),
        np.array(chg.chg).shape,
    ) == (
        False,
        False,
        (2, 56, 56, 56),
    )
    vac = loc.vac_potential()[0]
    # assert round(vac, 2) == round(7.62302803577618, 2)

    # td = loc.to_dict()
    # fd = Locpot.from_dict(td)

    vac = loc.vac_potential(direction="Y")[0]
    # assert round(vac, 2) == round(7.62302803577618, 2)

    vac = loc.vac_potential(direction="Z")[0]
    # assert round(vac, 2) == round(7.62302803577618, 2)


def test_vrun():
    # print ('gapp',round(vrun.get_indir_gap,2))
    assert (round(vrun.get_indir_gap[0], 2)) == (0.73)
    gap2 = vrun.bandgap_occupation_tol()
    assert round(gap2[0], 2) == 0.73
    assert (round(vrun.get_dir_gap, 2)) == (2.62)
    vrun.get_bandstructure(kpoints_file_path=band_kp, plot=True)
    assert (round(opt_vrun.get_dir_gap, 2)) == (2.62)
    assert (vrun.total_dos[0][0]) == -8.1917
    assert vrun.converged == True
    # TODO Serious issue: assert (opt_vrun.total_dos[0][0]) == -8.1917
    assert (vrun.eigenvalues[0][0][0][0]) == -6.1917
    assert (opt_vrun.eigenvalues[0][0][0][0]) == -6.1917
    assert vrun.is_spin_polarized == True
    assert opt_vrun.is_spin_polarized == False
    assert vrun.is_spin_orbit == False
    assert list(vrun.fermi_velocities[0]) == []
    pdos1 = vrun.partial_dos_spdf
    pdos2 = vrun.projected_atoms_spins_kpoints_bands
    pdos3 = vrun.projected_spins_kpoints_bands
    pdos4 = vrun.get_atom_resolved_dos(plot=True)
    pdos5 = vrun.get_spdf_dos(plot=True)
    pdos1 = opt_vrun.partial_dos_spdf
    pdos2 = opt_vrun.projected_atoms_spins_kpoints_bands
    pdos3 = opt_vrun.projected_spins_kpoints_bands
    # TODO pdos4 = opt_vrun.get_atom_resolved_dos()
    # TODO #pdos5 = opt_vrun.get_spdf_dos()
    td = vrun.to_dict()
    fd = Vasprun.from_dict(td)
    vrun_dm = Vasprun(
        os.path.join(os.path.dirname(__file__), "JVASP-86924.xml")
    )
    fv = vrun_dm.fermi_velocities
    assert round(fv[0][0], 2) == round(491630.23058338434, 2)
    # Based on bug fixes
    tdos = non_spinpol_vrun.total_dos
    diel_op_x, diel_op_y = vasp544_mbj_optics_vrun.dielectric_loptics


def test_osz():
    assert (float(osz.magnetic_moment)) == (0.0)
    assert osz.electronic_steps[0][2] == "0.747368666078E+01"
    td = osz.to_dict()
    fd = Oszicar.from_dict(td)


def test_out():
    assert (round(out.elastic_props()["KV"], 2)) == (87.27)
    out_efg = Outcar(
        os.path.join(os.path.dirname(__file__), "OUTCAR.EFG-JVASP-12148")
    )
    print(out_efg.all_structures())
    rl, imag = opt_out.freq_dielectric_tensor()
    nedos = opt_out.nedos
    out_efg_raw = Outcar(
        os.path.join(os.path.dirname(__file__), "OUTCAR.EFG-JVASP-12148")
    ).efg_raw_tensor
    assert out_efg.efg_tensor_diag()[0][0] == -4.766
    assert out_efg.efg_tensor_diag(std_conv=False)[0][0] == -4.766
    assert out_efg.quad_mom[0][0] == 0.023
    assert out_efg.converged == True
    td = out_efg.to_dict()
    fd = Outcar.from_dict(td)
    print("out_efg_raw", (out_efg_raw))
    print()
    print("out_efg_raw", np.linalg.eig(out_efg_raw)[0])
    print()
    print("out_efg", out_efg.efg_tensor_diag())
    print()
    # TODO: compare withvasprun gap
    gap = out_efg.bandgap


def test_dfpt():
    vrun = Vasprun(
        os.path.join(os.path.dirname(__file__), "vasprun.xml.JVASP-39")
    )
    out = Outcar(os.path.join(os.path.dirname(__file__), "OUTCAR.JVASP-39"))
    bec = round(vrun.dfpt_data["born_charges"][0][0][0], 2)
    eig = round(out.phonon_eigenvalues[2], 2)
    ionic_pz, total_pz = out.piezoelectric_tensor
    pz = total_pz[2][0]
    assert (bec, eig, pz) == (2.52, 19.58, -0.26756)
    # print (vrun.all_stresses)
    assert vrun.all_stresses[0][0][0] == -14.79381147
    assert vrun.all_forces[0][0][0] == 0
    assert vrun.all_energies[0] == -24.86360178
    assert vrun.all_structures[0].volume == 42.60334334259966


def single_element_vrun():
    vrun = Vasprun(
        os.path.join(os.path.dirname(__file__), "vasprun.xml.JVASP-816")
    )
    p = vrun.partial_dos_spdf()


# def test_waveder():
#    assert (
#        np.iscomplex(wder.get_orbital_derivative_between_states(0, 0, 0, 0, 0))
#        == True
#    )
#    assert (
#        complex(wder.get_orbital_derivative_between_states(0, 0, 0, 0, 0)).real
#    ) == -2.216161544844864e-15
#    assert (wder.nbands, wder.nkpoints, wder.nelect) == (36, 56, 8)


def test_ir():
    vrun = Vasprun(
        os.path.join(os.path.dirname(__file__), "vasprun.xml.JVASP-39")
    )
    out = Outcar(os.path.join(os.path.dirname(__file__), "OUTCAR.JVASP-39"))
    phonon_eigenvectors = vrun.dfpt_data["phonon_eigenvectors"]
    vrun_eigs = vrun.dfpt_data["phonon_eigenvalues"]
    phonon_eigenvalues = out.phonon_eigenvalues
    masses = vrun.dfpt_data["masses"]
    born_charges = vrun.dfpt_data["born_charges"]
    x, y = ir_intensity(
        phonon_eigenvectors=phonon_eigenvectors,
        phonon_eigenvalues=phonon_eigenvalues,
        masses=masses,
        born_charges=born_charges,
    )
    assert round(y[2], 2) == 0


def test_wavecar():
    gvec = wf_noso.gvectors()
    assert (gvec.shape) == (555, 3)


def test_raman():
    ram = os.path.join(os.path.dirname(__file__), "vasp_raman.dat")
    parse_raman_dat(ram)


# test_out()
