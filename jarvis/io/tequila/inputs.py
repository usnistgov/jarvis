"""Methods for handling input/output files for tequilla."""
# Reference: https://arxiv.org/abs/2102.11452
import numpy
import numpy as np
import tequila as tq
from jarvis.io.qiskit.inputs import HermitianSolver


def get_eigvals(mat=[]):
    """Get eigenvalues for Hermitianmatrix."""
    herm = HermitianSolver(mat)
    matrix = herm.mat
    H = tq.paulis.Zero()
    rows = matrix.shape[0]
    columns = matrix.shape[0]
    n_qubits = herm.n_qubits()
    print("true values", np.linalg.eigh(matrix)[0])
    for i in range(rows):
        for j in range(columns):
            H += matrix[i, j] * tq.paulis.KetBra(
                ket=i, bra=j, n_qubits=n_qubits
            )

    hermitian, anti = H.split()
    eigenValues, eigenVectors = numpy.linalg.eigh(hermitian.to_matrix())
    # print(tq.QubitWaveFunction(eigenVectors[:, 2]))

    circuits = []
    energies = []
    factor = 22
    # TODO: replace factor with reverse VQE value
    opt_variables = {}
    for i in range(rows):
        U = tq.gates.Ry(angle=(0, i), target=0)
        U += tq.gates.Ry(angle=(1, i), target=1)
        U += tq.gates.Ry(angle=(2, i), target=2)
        active_variables = U.extract_variables()
        E = tq.ExpectationValue(H=hermitian, U=U)
        ex_objective = E
        P1 = tq.paulis.Projector("1.0*|000>")
        for k in range(i):
            S = tq.ExpectationValue(H=P1, U=U + circuits[k].dagger())
            ex_objective += factor * abs(energies[k]) * S
        opt_variables = {
            **opt_variables,
            **{k: 1.0e-3 for k in active_variables},
        }
        result = tq.minimize(
            objective=ex_objective,
            method="bfgs",
            variables=active_variables,
            initial_values=opt_variables,
        )
        circuits.append(U)
        energies.append(result.energy)
        opt_variables = {**opt_variables, **result.variables}
    return energies


# def su2_circuit(n_qubits, layers, label):
#     """Get  SU2 circuit."""
#     U = tq.QCircuit()
#     for layer in range(layers):
#         for n in range(n_qubits):
#             U += tq.gates.Rz(target=n, angle=("z", n, layer, label))
#             U += tq.gates.Ry(target=n, angle=("y", n, layer, label))
#             If n < n_qubits-1:
#                 U += tq.gates.X(control=n, target=n+1)
#     return U

"""
if __name__ == "__main__":
    from jarvis.db.figshare import (
        get_wann_phonon,
        get_hk_tb,
        get_wann_electron,
    )

    w, atoms, o = get_wann_electron("JVASP-816")
    # w, atoms = get_wann_phonon("JVASP-816")
    hk = get_hk_tb(w=w, k=[0.0, 0.0, 0])
    print("hk", hk)
    H = HermitianSolver(hk)
    ens = get_eigvals(hk)
    print(ens)
"""
