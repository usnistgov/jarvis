from jarvis.spillage.spillage import *
import os

wf_so = os.path.join(os.path.dirname(__file__), "..", "spillage", "WAVECAR.soc")

wf_noso = os.path.join(os.path.dirname(__file__), "..", "spillage", "WAVECAR.nosoc")


def test_spillage():
    gamma_max, gamma_k, kmax, kpoints, noso_direct, so_direct, x, y, so_lumo, so_homo, noso_lumo, noso_homo = overlap_so_spinpol(
        wf_noso, wf_so
    )
    assert gamma_max == 1.3634111271008775
