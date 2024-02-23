"""Module to solve Hermitian Matrix and predict badstructures."""
# Reference: https://doi.org/10.1088/1361-648X/ac1154

import numpy as np
import itertools
import functools
from qiskit_aer import Aer
from qiskit.utils import QuantumInstance, algorithm_globals
from qiskit.opflow import I, X, Y, Z
from qiskit.algorithms import VQE
from jarvis.core.kpoints import Kpoints3D as Kpoints
from jarvis.db.figshare import get_hk_tb
from jarvis.core.kpoints import generate_kgrid
import matplotlib.pyplot as plt

plt.switch_backend("agg")
# from qiskit.algorithms.optimizers import SLSQP


def decompose_Hamiltonian(H):
    """Decompose Hermitian matrix into Pauli basis."""
    # Inspired from
    # https://github.com/PennyLaneAI/pennylane/blob/master/pennylane/utils.py#L45
    # https://qiskit.org/documentation/tutorials/algorithms/04_vqe_advanced.html
    x, y = H.shape
    N = int(np.log2(len(H)))
    if len(H) - 2 ** N != 0 or x != y:
        raise ValueError(
            "Hamiltonian should be in the form (2^n x 2^n), for any n>=1"
        )
    pauilis = [I, X, Y, Z]
    decomposedH = 0
    for term in itertools.product(pauilis, repeat=N):
        matrices = [i.to_matrix() for i in term]
        # coefficient of the pauli string = (1/2^N) * (Tr[pauliOp x H])
        coeff = np.trace(functools.reduce(np.kron, matrices) @ H) / (2 ** N)
        coeff = np.real_if_close(coeff).item()
        if coeff == 0:
            continue
        obs = 1
        for i in term:
            obs = obs ^ i
        decomposedH += coeff * obs
    return decomposedH


class HermitianSolver(object):
    """Solve a Hermitian matrix using quantum algorithms."""

    def __init__(self, mat=[], verbose=False):
        """Initialize with a numpy Hermitian matrix."""
        N = int(np.ceil(np.log2(len(mat))))
        hk = np.zeros((2 ** N, 2 ** N), dtype="complex")
        hk[: mat.shape[0], : mat.shape[1]] = mat
        self.mat = hk
        self.verbose = verbose
        if not self.check_hermitian():
            raise ValueError("Only implemented for Hermitian matrix.")

    def n_qubits(self):
        """Get number of qubits required."""
        return int(np.log2(len(self.mat)))

    def check_hermitian(self):
        """Check if a matrix is Hermitian."""
        adjoint = self.mat.conj().T
        return np.allclose(self.mat, adjoint)

    def run_vqe(
        self,
        backend=Aer.get_backend("statevector_simulator"),
        var_form=None,
        optimizer=None,
        reps=None,
        mode="min_val",
    ):
        """Run variational quantum eigensolver."""
        seed = 50
        algorithm_globals.random_seed = seed
        N = self.n_qubits()
        qi = QuantumInstance(
            Aer.get_backend("statevector_simulator"),
            seed_transpiler=seed,
            seed_simulator=seed,
        )
        # n_qubits = self.n_qubits

        if mode == "max_val":
            Hamil_qop = decompose_Hamiltonian(-1 * self.mat)
            np_eig = min(np.linalg.eig(-1 * self.mat)[0])
            if self.verbose:
                print("np_eig", np_eig)
        else:
            Hamil_qop = decompose_Hamiltonian(self.mat)
            np_eig = min(np.linalg.eig(self.mat)[0])
            if self.verbose:
                print("np_eig", np_eig)

        if var_form is None:
            if reps is None:
                reps = 2
                # reps=5
            from qiskit.circuit.library import EfficientSU2

            var_form = EfficientSU2(N, reps=reps)
        if optimizer is None:
            vqe = VQE(var_form, quantum_instance=qi)

        else:
            vqe = VQE(var_form, optimizer=optimizer, quantum_instance=qi)
        result = vqe.compute_minimum_eigenvalue(operator=Hamil_qop)
        en = result.eigenvalue

        if mode == "max_val":
            en = -1 * en

        return en, result, vqe

    def run_numpy(self):
        """Obtain eigenvalues and vecs using Numpy solvers."""
        return np.linalg.eigh(self.mat)

    # def run_qpe(self, n_ancillae=8):
    #    """Run quantum phase estimations."""
    #    quantum_instance = aqua.QuantumInstance(
    #        backend=Aer.get_backend("statevector_simulator"), shots=1
    #    )
    #    Hamil_mat = aqua.operators.MatrixOperator(self.mat)
    #    # Hamil_mat = MatrixOperator(self.mat)
    #    #Hamil_qop = aqua.operators.op_converter.to_weighted_pauli_operator(
    #    #    Hamil_mat
    #    #)
    #    # Hamil_qop = op_converter.to_weighted_pauli_operator(Hamil_mat)

    #    Hamil_qop = decompose_Hamiltonian(self.mat)
    #    qpe = aqua.algorithms.QPE(Hamil_qop, num_ancillae=n_ancillae)
    #    qpe_result = qpe.run(quantum_instance)
    #    # qc = qpe.construct_circuit(measurement=True)
    #    print("qpe_result", qpe_result)
    #    return qpe_result["eigenvalue"], qpe_result, qpe

    def run_vqd(
        self,
        backend=Aer.get_backend("statevector_simulator"),
        var_form=None,
        optimizer=None,
        reps=2,
        # reps=5,
    ):
        """Run variational quantum deflation."""
        tmp = HermitianSolver(self.mat)
        max_eigval, vqe_result, vqe = tmp.run_vqe(
            backend=backend,
            var_form=var_form,
            optimizer=optimizer,
            reps=reps,
            mode="max_val",
        )
        eigvals = [max_eigval]
        eigstates = [vqe_result.eigenstate]
        #         eigvals = []
        #         eigstates= []
        for r in range(len(tmp.mat) - 1):
            val, vqe_result, vqe = tmp.run_vqe(
                backend=backend,
                var_form=var_form,
                optimizer=optimizer,
                reps=reps,
            )
            outer_prod = np.outer(
                vqe_result.eigenstate, np.conj(vqe_result.eigenstate).T
            )
            tmp.mat = tmp.mat - (val - max_eigval) * outer_prod
            eigvals.append(val)
            eigstates.append(vqe_result.eigenstate)
            tmp = HermitianSolver(tmp.mat)

        eigvals = np.array(eigvals)
        eigstates = np.array(eigstates)
        order = np.argsort(eigvals)
        eigvals = eigvals[order]
        eigstates = eigstates[order]
        return eigvals, eigstates


def get_bandstruct(
    w=[],
    atoms={},
    ef=0,
    line_density=1,
    ylabel="eV",  # "Energy ($cm^{-1}$)",
    font=22,
    var_form=None,
    filename="bands.png",
    savefig=True,
    neigs=None,
    max_nk=None,
    tol=None,
    factor=1,
    verbose=False,
):
    """Compare bandstructures using quantum algos."""
    info = {}
    kpoints = Kpoints().kpath(atoms, line_density=line_density)
    labels = kpoints.to_dict()["labels"]
    kpts = kpoints.to_dict()["kpoints"]
    if verbose:
        print("Number of kpoints:", len(kpts))

    eigvals_q = []
    eigvals_np = []
    for ii, i in enumerate(kpts):
        if max_nk is not None and ii == max_nk:
            break
            # For reducing CI/CD time
            print("breaking here", ii, max_nk)
        else:
            try:

                hk = get_hk_tb(w=w, k=i)
                HS = HermitianSolver(hk)
                vqe_vals, _ = HS.run_vqd(var_form=var_form)
                np_vals, _ = HS.run_numpy()
                if verbose:
                    print("kp=", ii, i)
                    print("np_vals", np_vals)
                    print("vqe_vals", vqe_vals)
                eigvals_q.append(vqe_vals)
                eigvals_np.append(np_vals)
                # break
                if (
                    neigs is not None
                    and isinstance(neigs, int)
                    and neigs == len(eigvals_q)
                ):
                    break
            except Exception as exp:
                print(exp)
                pass
    eigvals_q = factor * np.array(eigvals_q)  # 3.14 for phonon
    eigvals_np = factor * np.array(eigvals_np)

    for ii, i in enumerate(eigvals_q.T - ef):
        if ii == 0:
            plt.plot(i, "*", c="b", label="VQD")
        else:
            plt.plot(i, "*", c="b")

    for ii, i in enumerate(eigvals_np.T - ef):
        if ii == 0:
            plt.plot(i, c="g", label="Numpy")
        else:
            plt.plot(i, c="g")
    new_kp = []
    new_labels = []
    count = 0
    kp = np.arange(len(kpts))
    for i, j in zip(kp, labels):
        if j != "":
            if count > 1 and count < len(labels) - 1:
                if labels[count] != labels[count + 1]:
                    new_kp.append(i)
                    new_labels.append("$" + str(j) + "$")
            else:
                new_kp.append(i)
                new_labels.append("$" + str(j) + "$")
        count += 1
    info["eigvals_q"] = list(eigvals_q.tolist())
    info["eigvals_np"] = list(eigvals_np.tolist())
    info["kpts"] = list(kpts)
    info["new_kp"] = list(np.array(new_kp).tolist())
    info["new_labels"] = list(new_labels)
    info["ef"] = ef
    if verbose:
        print(info)
    if tol is not None:
        plt.ylim([tol, np.max(eigvals_q)])
    plt.rcParams.update({"font.size": font})
    plt.xticks(new_kp, new_labels)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()

    if savefig:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()
    return info


def get_dos(
    w=[],
    grid=[2, 1, 1],
    proj=None,
    efermi=0.0,
    xrange=None,
    nenergy=100,
    sig=0.02,
    use_dask=True,
    filename="dos.png",
    savefig=True,
):
    """Get density of states."""
    nwan = int(np.ceil(np.log2(w.nwan))) ** 2
    kpoints = generate_kgrid(grid=grid)
    nk = len(kpoints)
    q_vals = np.zeros((nk, nwan), dtype=float)
    np_vals = np.zeros((nk, nwan), dtype=float)
    pvals = np.zeros((nk, nwan - 1), dtype=float)
    # if use_dask:
    # def get_vqd_vals(k):
    #    hk = get_hk_tb(w=w, k=k)
    #    HS = HermitianSolver(hk)
    #    vqe_vals, _ = HS.run_vqd()
    #    return vqe_vals

    # values=[delayed(get_vqd_vals)(k) for k in  kpoints]
    # resultsDask = compute(*values, scheduler='processes')
    for i, k in enumerate(kpoints):
        hk = get_hk_tb(w=w, k=k)
        HS = HermitianSolver(hk)
        vqe_vals, _ = HS.run_vqd()
        n_vals, _ = HS.run_numpy()
        print("np_vals", n_vals, len(n_vals), np_vals.shape)
        print("vqe_vals", vqe_vals, len(vqe_vals), q_vals.shape)
        q_vals[i, :] = vqe_vals
        np_vals[i, :] = n_vals

    if xrange is None:
        vmin = np.min(q_vals[:])
        vmax = np.max(q_vals[:])
        vmin2 = vmin - (vmax - vmin) * 0.05
        vmax2 = vmax + (vmax - vmin) * 0.05
        xrange = [vmin2, vmax2]
        # plt.xlim(xrange)

    energies = np.arange(
        xrange[0], xrange[1] + 1e-5, (xrange[1] - xrange[0]) / float(nenergy),
    )
    dos = np.zeros(np.size(energies))
    pdos = np.zeros(np.size(energies))

    v = q_vals

    #   condmin = np.min(v[v > 0.0])
    #   valmax = np.max(v[v < 0.0])
    #   print("DOS BAND GAP ", condmin - valmax, "    ", valmax, " ", condmin)

    c = -0.5 / sig ** 2
    for i in range(np.size(energies)):
        arg = c * (v - energies[i]) ** 2
        dos[i] = np.sum(np.exp(arg))
        if proj is not None:
            pdos[i] = np.sum(np.exp(arg) * pvals)

    de = energies[1] - energies[0]
    dos = dos / sig / (2.0 * np.pi) ** 0.5 / float(nk)
    if proj is not None:
        pdos = pdos / sig / (2.0 * np.pi) ** 0.5 / float(nk)
    print("np.sum(dos) ", np.sum(dos * de))
    if proj is not None:
        print("np.sum(pdos) ", np.sum(pdos * de))
    plt.plot(energies, dos)
    if savefig:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()
    return energies, dos, pdos


"""
if __name__ == "__main__":
    from jarvis.db.figshare import (
        get_wann_electron,
        get_wann_phonon,
    )
    from jarvis.core.circuits import QuantumCircuitLibrary

    wtbh, ef, atoms = get_wann_electron("JVASP-816")
    hk = get_hk_tb(w=wtbh, k=[0.0, 0.0, 0.0])
    H = HermitianSolver(hk)
    qc = QuantumCircuitLibrary(n_qubits=3).circuit6()  # 2^3  = 8
    en, vqe_result, vqe = H.run_vqe(mode="min_val", var_form=qc)
"""
