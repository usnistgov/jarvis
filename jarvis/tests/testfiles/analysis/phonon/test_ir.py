from jarvis.io.vasp.outputs import Vasprun, Outcar
from jarvis.analysis.phonon.ir import ir_intensity

import os

out = Outcar(
    os.path.join(
        os.path.dirname(__file__), "..", "..", "io", "vasp", "OUTCAR.JVASP-39"
    )
)
vrun = Vasprun(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "io",
        "vasp",
        "vasprun.xml.JVASP-39",
    )
)


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
    assert max(y) == 0.3511482090386446
    pdos = vrun.partial_dos_spdf


test_ir()
