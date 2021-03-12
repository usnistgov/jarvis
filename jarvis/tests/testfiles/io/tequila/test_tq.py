"""
# Excluding tests because tequila doesn't seem to install well.
from jarvis.io.tequila.inputs import get_eigvals
from jarvis.io.qiskit.inputs import HermitianSolver
from jarvis.db.figshare import (
    get_wann_phonon,
    get_hk_tb,
    get_wann_electron,
)

def test_inp():
    w, atoms, o = get_wann_electron("JVASP-816")
    # w, atoms = get_wann_phonon("JVASP-816")
    hk = get_hk_tb(w=w, k=[0.0, 0.0, 0])
    print("hk", hk)
    H = HermitianSolver(hk)
    ens = get_eigvals(hk)
    print(ens)

"""
