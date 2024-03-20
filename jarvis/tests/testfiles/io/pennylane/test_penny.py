"""
# Commenting due to pypi conflict in qiskit
from jarvis.io.pennylane.inputs import run_vqe
from jarvis.db.figshare import (
    get_wann_phonon,
    get_hk_tb,
    get_wann_electron,
)

def test_inp():
    w, atoms, o = get_wann_electron("JVASP-816")
    # w, atoms = get_wann_phonon("JVASP-816")
    hk = get_hk_tb(w=w, k=[0.0, 0.0, 0])
    run_vqe(mat=hk)
"""
