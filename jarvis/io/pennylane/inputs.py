"""Methods for handling input/output files for pennylane."""
# Reference: https://arxiv.org/abs/2102.11452
from pennylane.utils import decompose_hamiltonian
import pennylane as qml
import numpy as np
from jarvis.io.qiskit.inputs import HermitianSolver


def run_vqe(mat=[], max_iterations=100, conv_tol=1e-04, step_size=0.01):
    """Run pennylane VQE solver."""
    herm = HermitianSolver(mat)
    hk = herm.mat
    nwires = herm.n_qubits()
    print("nwires=", nwires)
    coeff, obs_list = decompose_hamiltonian(hk)
    H = qml.Hamiltonian(coeff, obs_list)
    # print(H)
    ansatz = qml.templates.StronglyEntanglingLayers
    dev = qml.device("default.qubit", wires=nwires)
    qml.enable_tape()
    cost_fn = qml.ExpvalCost(ansatz, H, dev)
    opt = qml.GradientDescentOptimizer(stepsize=step_size)
    params = np.random.rand(
        1, nwires, 3
    )  # torch.rand([2, nwires, 3]).cpu().detach().numpy() #init_params
    prev_energy = cost_fn(params)
    gd_param_history = [params]
    gd_cost_history = [prev_energy]
    for n in range(max_iterations):
        # Take step
        params = opt.step(cost_fn, params)
        gd_param_history.append(params)
        # print (n,ansatz)
        # print ()
        # Compute energy
        energy = cost_fn(params)
        gd_cost_history.append(energy)
        # Calculate difference between new and old energies
        conv = np.abs(energy - prev_energy)
        if n % 20 == 0:
            print(
                "Iteration={:}, Energy={:.8f} Ha,Convergence parameter = {"
                ":.8f} Ha".format(n, energy, conv)
            )
        if conv <= conv_tol:
            break
        prev_energy = energy

    # print()
    print("Final value of the energy = {:.8f} Ha".format(energy))
    # print("Number of iterations = ", n)
    vect = np.array(cost_fn.hamiltonian.coeffs).reshape(hk.shape)
    return energy, vect


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
    run_vqe(mat=hk)
"""
