from jarvis.io.vasp.outputs import Vasprun, Outcar
from jarvis.analysis.phonon.ir import ir_intensity, ir_intensity_phonopy
import os
import shutil

out = Outcar(
    os.path.join(
        os.path.dirname(__file__), "..", "..", "io", "vasp", "OUTCAR.JVASP-39"
    )
)
vrun_file = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "io",
    "vasp",
    "vasprun.xml.JVASP-39",
)
dirc = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "io",
    "vasp",
)

vrun = Vasprun(vrun_file)


def test_ir():
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
    print(max(y))
    assert round(max(y), 2) == round(0.3511482090386446, 2)
    pdos = vrun.partial_dos_spdf
    cwd = os.getcwd()
    if not os.path.exists("vasprun.xml"):
        shutil.copy2(vrun_file, "vasprun.xml")
    x, y = ir_intensity_phonopy(vasprun=vrun_file)  # , run_dir=dirc)


test_ir()
