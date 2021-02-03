#
#from qiskit.circuit.library import EfficientSU2
from qiskit.circuit import QuantumCircuit
from qiskit.circuit import library
from qiskit import circuit
from jarvis.db.figshare import (
    get_wann_phonon,
    get_hk_tb,
    get_wann_electron,
)
from jarvis.core.atoms import Atoms
from jarvis.db.jsonutils import dumpjson
from jarvis.io.qiskit.inputs import HermitianSolver,get_bandstruct,get_dos
import matplotlib.pyplot as plt
plt.switch_backend('agg')

def test_inp():
    w, ef, atoms = get_wann_electron("JVASP-816")
#     info = get_bandstruct(
#         w=w,
#         line_density=1,
#         atoms=atoms,
#         ef=ef,
#         filename="Alelect.png",
#         ylabel="Energy (eV)",
#     )
    #dumpjson(data=info, filename="Alelect.json")
    w, atoms = get_wann_phonon("JVASP-816", factor=34.3)
    hk = get_hk_tb(w=w, k=[0.0, 0.0, 0.0])
    var_form = QuantumCircuit(2)
    H = HermitianSolver(hk)
    en, vqe_result, vqe = H.run_vqe(mode="max_val", var_form=var_form,optimizer=None)
    print("en=", en)
    #info = get_bandstruct(
    #    w=w,
    #    line_density=1,
    #    atoms=atoms,
    #    tol=0.1,
    #    filename="Alphon.png",
    #    ylabel="Freq.(cm$^{-1}$)",
    #)
    ## dumpjson(data=info,filename='Alphon.json')
    #eigs,vecs=H.run_vqd()
    # print(eigs)
    # print(vecs)

    eigs, vecs = H.run_numpy()
    print(eigs)
    # print(vecs)
    # get_bandstruct(w=w, atoms=atoms, tol=0.1)
    #get_dos(w=w)
    #H.run_qpe()

