"""Module to design quantum circuits with qiskit."""
from qiskit.circuit import QuantumCircuit, ParameterVector


def tb_ansatz(n_qubits=9, reps=2):
    """Generate tight-binding ansatz."""
    circ = QuantumCircuit(n_qubits)
    params = ParameterVector("theta", reps * 4 * n_qubits)
    count = 0
    for r in range(reps):

        for i in range(n_qubits):
            circ.ry(params[count], i)
            count += 1
            circ.rz(params[count], i)
            count += 1
        for i in range(n_qubits):

            if i != n_qubits - 1:
                circ.cx(i, n_qubits - 1)

            for j in range(n_qubits):
                if i + j < n_qubits and j > 0:
                    circ.cx(i, i + j)
        for i in range(n_qubits):
            circ.ry(params[count], i)
            count += 1
            circ.rz(params[count], i)
            count += 1
    return circ
