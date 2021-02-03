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
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit import aqua
import numpy as np
params = np.random.rand(3)
optimizer = aqua.components.optimizers.COBYLA(maxiter=500, tol=0.0001)
params = np.random.rand(3)
from qiskit import Aer, execute
backend = Aer.get_backend("qasm_simulator")
NUM_SHOTS = 10000
# Obtain the output distribution using the final parameters

target_distr = np.random.rand(2)
# We now convert the random vector into a valid probability vector
target_distr /= sum(target_distr)

def get_var_form(params):
    qr = QuantumRegister(3, name="q")
    cr = ClassicalRegister(1, name='c')
    qc = QuantumCircuit(qr, cr)
    qc.u2(params[0], params[1], qr[0])#, qr[0])
    qc.measure(qr, cr[0])
    return qc
def get_probability_distribution(counts):
    output_distr = [v / NUM_SHOTS for v in counts.values()]
    if len(output_distr) == 1:
        output_distr.append(1 - output_distr[0])
    return output_distr

def objective_function(params):
    # Obtain a quantum circuit instance from the paramters
    qc = get_var_form(params)
    # Execute the quantum circuit to obtain the probability distribution associated with the current parameters
    result = execute(qc, backend, shots=NUM_SHOTS).result()
    # Obtain the counts for each measured state, and convert those counts into a probability vector
    output_distr = get_probability_distribution(result.get_counts(qc))
    # Calculate the cost as the distance between the output distribution and the target distribution
    cost = sum([np.abs(output_distr[i] - target_distr[i]) for i in range(2)])
    return cost

ret = optimizer.optimize(num_vars=3, objective_function=objective_function, initial_point=params)
qc = get_var_form(ret[0])

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
    var_form = qc
    #var_form.measure([0, 1], [0, 1])
    H = HermitianSolver(hk)
    en, vqe_result, vqe = H.run_vqe(mode="max_val", var_form=var_form,optimizer=optimizer)
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

