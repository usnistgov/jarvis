from jarvis.io.vasp.outputs import Vasprun, Oszicar, Wavecar, Waveder, Chgcar, Outcar
import numpy as np
import os
from jarvis.analysis.phonon.ir import ir_intensity
import matplotlib.pyplot as plt
plt.switch_backend('agg')
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
        "MAIN-BAND-bulk@mp_149",
        "vasprun.xml",
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
        "MAIN-BAND-bulk@mp_149",
        "KPOINTS",
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
        "MAIN-ELASTIC-bulk@mp_149",
        "OUTCAR",
    )
)
wder = Waveder(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "MAIN-OPTICS-bulk@mp_149",
        "WAVEDER",
    )
)
wf_noso = Wavecar(
    filename=os.path.join(
        os.path.dirname(__file__), "..", "..", "analysis", "spillage", "WAVECAR.nosoc"
    ),
    lsorbit=False,
)
wf_so = Wavecar(
    filename=os.path.join(
        os.path.dirname(__file__), "..", "..", "analysis", "spillage", "WAVECAR.soc"
    ),
    lsorbit=True,
)


def test_chgcar():
    # print (chg.is_spin_polarized(), chg.is_spin_orbit(), np.array(chg.chg).shape)
    assert (chg.is_spin_polarized(), chg.is_spin_orbit(), np.array(chg.chg).shape) == (
        True,
        False,
        (2, 56, 56, 56),
    )


def test_vrun():
    # print ('gapp',round(vrun.get_indir_gap,2))
    assert (round(vrun.get_indir_gap, 2)) == (0.73)
    assert (round(vrun.get_dir_gap, 2)) == (2.62)
    vrun.get_bandstructure(kpoints_file_path = band_kp)
    assert (round(opt_vrun.get_dir_gap, 2)) == (2.62)
    assert (vrun.total_dos[0][0]) == -8.1917
    # TODO Serious issue: assert (opt_vrun.total_dos[0][0]) == -8.1917

def test_osz():
    assert (float(osz.magnetic_moment)) == (0.0)


def test_out():
    assert (round(out.elastic_props()["KV"], 2)) == (87.27)
    out_efg = Outcar(os.path.join(os.path.dirname(__file__), "OUTCAR.EFG-JVASP-12148"))
    assert out_efg.efg_tensor_diag[0][0] == -4.766
    assert out_efg.quad_mom[0][0] == 0.023
    assert out_efg.converged == True

def test_dfpt():
   vrun = Vasprun(os.path.join(os.path.dirname(__file__), "vasprun.xml.JVASP-39"))
   out = Outcar(os.path.join(os.path.dirname(__file__), "OUTCAR.JVASP-39"))
   bec = round(vrun.dfpt_data['born_charges'][0][0][0],2)
   eig = round(out.phonon_eigenvalues[2],2)
   ionic_pz,total_pz = out.piezoelectric_tensor
   pz = total_pz[2][0]
   assert (bec, eig, pz)==(2.52, 19.58,-0.26756) 
   # print (vrun.all_stresses)
   assert vrun.all_stresses[0][0][0]==-14.79381147
def test_waveder():
   assert np.iscomplex(wder.get_orbital_derivative_between_states(0,0,0,0,0)) == True
   assert (complex(wder.get_orbital_derivative_between_states(0,0,0,0,0)).real) == -2.216161544844864e-15
   assert (wder.nbands, wder.nkpoints, wder.nelect) == (36, 56, 8)
def test_ir():
    vrun = Vasprun(os.path.join(os.path.dirname(__file__), "vasprun.xml.JVASP-39"))
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
    assert round(y[2],2)==1.38

def test_wavecar():
   gvec = wf_noso.gvectors() 
   assert (gvec.shape)==(555,3)


