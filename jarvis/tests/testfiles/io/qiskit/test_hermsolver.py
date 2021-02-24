#
#from qiskit.circuit.library import EfficientSU2
from qiskit.circuit import QuantumCircuit
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
from qiskit.circuit import QuantumCircuit, ParameterVector
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit import aqua
import numpy as np
from qiskit import Aer, execute

import math

def nCr(n,r):
    f = math.factorial
    return int(f(n) / f(r) / f(n-r))

def variational_circuit(num_qubits = 2,reps = 1):
    # import required qiskit libraries if additional libraries are required
    # build the variational circuit
    #var_circuit = EfficientSU2(num_qubits=3, su2_gates= ['rx', 'ry'], entanglement='circular', reps=3)
    #var_circuit = EfficientSU2(num_qubits=4, su2_gates= ['rx', 'ry'], entanglement='circular', reps=3)
    
    # return the variational circuit which is either a VaritionalForm or QuantumCircuit object
    from qiskit.circuit import QuantumCircuit, ParameterVector
    x = ParameterVector('x', length=num_qubits)  # creating a list of Parameters
    custom_circ = QuantumCircuit(num_qubits)

    # defining our parametric form
    for _ in range(reps):
        for i in range(num_qubits):
            custom_circ.rx(x[i], i)
        for i in range(num_qubits):
            for j in range(i + 1, num_qubits):
                custom_circ.cx(i, j)
                custom_circ.u1(x[i] * x[j], j)
                custom_circ.cx(i, j)
            
            
    
                
            

qc=variational_circuit()

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
    print ('shape=',hk.shape)
    H = HermitianSolver(hk)
    en, vqe_result, vqe = H.run_vqe(mode="max_val", var_form=qc)#,optimizer=optimizer)
    print("en=", en)
    info = get_bandstruct(
        w=w,
        line_density=1,
        atoms=atoms,
        tol=0.1,
        neigs=2,
        max_nk=2,
        filename="Alphon.png",
        ylabel="Freq.(cm$^{-1}$)",
    )
    ## dumpjson(data=info,filename='Alphon.json')
    #eigs,vecs=H.run_vqd()
    # print(eigs)
    # print(vecs)

    eigs, vecs = H.run_numpy()
    print(eigs)
    # print(vecs)
    # get_bandstruct(w=w, atoms=atoms, tol=0.1)
    get_dos(w=w)
    H.run_qpe()
