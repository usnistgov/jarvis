"""Module to design quantum circuits with qiskit."""
from qiskit.circuit import QuantumCircuit, ParameterVector

from qiskit.circuit.library import EfficientSU2, PauliTwoDesign, RealAmplitudes


class QuantumCircuitLibrary(object):
    """Module for storing various quantum circuits."""

    def __init__(self, n_qubits=3, reps=1):
        """Initialize class."""
        self.n_qubits = n_qubits
        self.reps = reps

    def circuit1(self):
        """Generate tight-binding ansatz."""
        reps = self.reps
        n_qubits = self.n_qubits
        circ = QuantumCircuit(n_qubits)
        params = ParameterVector("θ", reps * 1 * n_qubits)
        count = 0
        for r in range(reps):
            for i in range(n_qubits):
                circ.ry(params[count], i)
                count += 1
        return circ

    def circuit2(self):
        """Generate tight-binding ansatz."""
        reps = self.reps
        n_qubits = self.n_qubits
        circ = QuantumCircuit(n_qubits)
        params = ParameterVector("θ", reps * 2 * n_qubits)
        count = 0
        for r in range(reps):

            for i in range(n_qubits):
                circ.ry(params[count], i)
                count += 1
                circ.rz(params[count], i)
                count += 1
        return circ

    def circuit3(self):
        """Generate tight-binding ansatz."""
        reps = self.reps
        n_qubits = self.n_qubits
        circ = QuantumCircuit(n_qubits)
        params = ParameterVector("θ", reps * 2 * n_qubits)
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

        return circ

    def circuit4(self):
        """Generate RealAmplitudes ansatz."""
        reps = self.reps
        n_qubits = self.n_qubits
        return RealAmplitudes(num_qubits=n_qubits, reps=reps)

    def circuit5(self):
        """Generate PauliTwoDesign ansatz."""
        reps = self.reps
        n_qubits = self.n_qubits
        return PauliTwoDesign(num_qubits=n_qubits, reps=reps)

    def circuit6(self):
        """Generate EfficientSU2 ansatz."""
        reps = self.reps
        n_qubits = self.n_qubits
        return EfficientSU2(num_qubits=n_qubits, reps=reps)

    def circuit7(self):
        """Generate tight-binding ansatz."""
        reps = self.reps
        n_qubits = self.n_qubits
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


"""
if __name__ == "__main__":
    from jarvis.db.figshare import (
        get_wann_phonon,
        get_hk_tb,
        get_wann_electron,
    )
    from jarvis.io.qiskit.inputs import HermitianSolver

    w, atoms, ef = get_wann_electron("JVASP-816")
    kps = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    for kk, k in enumerate(kps):
        print()
        hk = get_hk_tb(w=w, k=k)
        HS = HermitianSolver(hk)
        n_qubits = HS.n_qubits()
        Q = QuantumCircuitLibrary(n_qubits=n_qubits)
        circs = [
            Q.circuit1(),
            Q.circuit2(),
            Q.circuit3(),
            Q.circuit4(),
            Q.circuit5(),
            Q.circuit6(),
        ]
        for ii, i in enumerate(circs):
            en, vqe_result, vqe = HS.run_vqe(var_form=i)
            vals, vecs = HS.run_numpy()
            diff = abs(en - vals[0])
            print("VQE,numpy", ii, kk, en, vals[0], diff)
            print(i)
"""
