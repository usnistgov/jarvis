from qiskit import Aer
from qiskit.utils import QuantumInstance, algorithm_globals
from qiskit.algorithms import VQE
from qiskit.algorithms.optimizers import SLSQP
import numpy as np
import itertools, functools
from qiskit.opflow import I, X, Y, Z
from jarvis.db.figshare import get_wann_electron, get_wann_phonon, get_hk_tb
from jarvis.core.circuits import QuantumCircuitLibrary
from jarvis.io.qiskit.inputs import HermitianSolver
from qiskit import Aer


def decompose_Hamiltonian(H):
    # Inspired from
    # https://github.com/PennyLaneAI/pennylane/blob/master/pennylane/utils.py#L45
    # https://qiskit.org/documentation/tutorials/algorithms/04_vqe_advanced.html
    x, y = H.shape
    N = int(np.log2(len(H)))
    if len(H) - 2**N != 0 or x != y:
        raise ValueError(
            "Hamiltonian should be in the form (2^n x 2^n), for any n>=1"
        )
    pauilis = [I, X, Y, Z]
    decomposedH = 0
    for term in itertools.product(pauilis, repeat=N):
        matrices = [i.to_matrix() for i in term]
        # coefficient of the pauli string = (1/2^N) * (Tr[pauliOp x H])
        coeff = np.trace(functools.reduce(np.kron, matrices) @ H) / (2**N)
        coeff = np.real_if_close(coeff).item()
        if coeff == 0:
            continue
        obs = 1
        for i in term:
            obs = obs ^ i
        decomposedH += coeff * obs
    return decomposedH


def test_qiskit():
    wtbh, Ef, atoms = get_wann_electron("JVASP-816")
    kpt = [0.5, 0.0, 0.5]  # X-point
    hk = get_hk_tb(w=wtbh, k=kpt)
    wtbh_op = decompose_Hamiltonian(hk)

    seed = 50
    algorithm_globals.random_seed = seed
    qi = QuantumInstance(
        Aer.get_backend("statevector_simulator"),
        seed_transpiler=seed,
        seed_simulator=seed,
    )
    n_qubits = int(np.log2(len(hk)))
    # ansatz = TwoLocal(rotation_blocks='ry', entanglement_blocks='cz')
    ansatz = QuantumCircuitLibrary(n_qubits=n_qubits, reps=1).circuit6()
    slsqp = SLSQP(maxiter=1000)
    vqe = VQE(ansatz, optimizer=slsqp, quantum_instance=qi)
    result = vqe.compute_minimum_eigenvalue(operator=wtbh_op)
    np_eig = min(np.linalg.eig(hk)[0])
    print("numpy min. eig", np_eig)

    eigenvalue = result.eigenvalue
    # print(result)
    print("VQE eig.", eigenvalue)


def test_statvector():
    backend = Aer.get_backend("statevector_simulator")
    # Aluminum JARVIS-ID: JVASP-816
    wtbh, Ef, atoms = get_wann_electron("JVASP-816")
    kpt = [0.5, 0.0, 0.5]  # X-point
    hk = get_hk_tb(w=wtbh, k=kpt)
    HS = HermitianSolver(hk)
    n_qubits = HS.n_qubits()
    circ = QuantumCircuitLibrary(n_qubits=n_qubits, reps=1).circuit6()
    en, vqe_result, vqe = HS.run_vqe(var_form=circ, backend=backend)
    vals, vecs = HS.run_numpy()
    # Ef: Fermi-level
    print("Classical, VQE (eV):", vals[0] - Ef, en - Ef)
    print("Show model\n", circ)


"""
# Commenting due to pypi conflicts in qiskit
#
# from qiskit.circuit.library import EfficientSU2
from qiskit.circuit import QuantumCircuit
from qiskit import circuit
from jarvis.db.figshare import (
    get_wann_phonon,
    get_hk_tb,
    get_wann_electron,
)
from jarvis.core.atoms import Atoms
from jarvis.db.jsonutils import dumpjson
from jarvis.io.qiskit.inputs import HermitianSolver, get_bandstruct, get_dos
import matplotlib.pyplot as plt

plt.switch_backend("agg")
from qiskit.circuit import QuantumCircuit, ParameterVector
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister

# from qiskit import aqua
import numpy as np
from qiskit import Aer, execute

import math
import gc

gc.collect()


def nCr(n, r):
    f = math.factorial
    return int(f(n) / f(r) / f(n - r))


# def test_hermitian_solver():
#     from jarvis.db.figshare import (
#         get_wann_electron,
#         get_wann_phonon,
#     )
#     from jarvis.core.circuits import QuantumCircuitLibrary

#     wtbh, ef, atoms = get_wann_electron("JVASP-816")
#     hk = get_hk_tb(w=wtbh, k=[0.0, 0.0, 0.0])
#     H = HermitianSolver(hk)
#     qc = QuantumCircuitLibrary(n_qubits=3).circuit6()  # 2^3  = 8
#     en, vqe_result, vqe = H.run_vqe(mode="min_val", var_form=qc)


def variational_circuit(num_qubits=2, reps=1):
    # import required qiskit libraries if additional libraries are required
    # build the variational circuit
    # var_circuit = EfficientSU2(num_qubits=3, su2_gates= ['rx', 'ry'], entanglement='circular', reps=3)
    # var_circuit = EfficientSU2(num_qubits=4, su2_gates= ['rx', 'ry'], entanglement='circular', reps=3)

    # return the variational circuit which is either a VaritionalForm or QuantumCircuit object
    from qiskit.circuit import QuantumCircuit, ParameterVector

    x = ParameterVector(
        "x", length=num_qubits
    )  # creating a list of Parameters
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


# qc = variational_circuit()


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
    # dumpjson(data=info, filename="Alelect.json")
    w, atoms = get_wann_phonon("JVASP-816", factor=34.3)
    hk = get_hk_tb(w=w, k=[0.0, 0.0, 0.0])
    print("shape=", hk.shape)
    H = HermitianSolver(hk)
    from jarvis.core.circuits import QuantumCircuitLibrary

    qc = QuantumCircuitLibrary(n_qubits=2).circuit1()
    # en, vqe_result, vqe = H.run_vqe(mode="max_val", reps=1)#,optimizer=optimizer)
    en, vqe_result, vqe = H.run_vqe(
        mode="max_val", var_form=qc
    )  # ,optimizer=optimizer)
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
    # eigs,vecs=H.run_vqd()
    # print(eigs)
    # print(vecs)

    eigs, vecs = H.run_numpy()
    print(eigs)
    # print(vecs)
    # get_bandstruct(w=w, atoms=atoms, tol=0.1)
    get_dos(w=w, grid=[2, 1, 1])
    # H.run_qpe()


#test_inp()
"""
