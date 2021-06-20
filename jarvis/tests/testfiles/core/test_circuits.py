from jarvis.core.circuits import QuantumCircuitLibrary


def test_circuits():
    Q = QuantumCircuitLibrary()
    circs = [
        Q.circuit1(),
        Q.circuit2(),
        Q.circuit3(),
        Q.circuit4(),
        Q.circuit5(),
        Q.circuit6(),
    ]
