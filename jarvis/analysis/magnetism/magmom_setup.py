"""Module for setting MAGMOM and AFM/FM orderings."""

from jarvis.analysis.structure.spacegroup import Spacegroup3D
import numpy as np


def check_match(a, b, tol=1e-4):
    """Check if a and b are the same, taking into account PBCs."""
    if abs(a[0] - b[0]) < tol or abs(abs(a[0] - b[0]) - 1) < tol:
        if abs(a[1] - b[1]) < tol or abs(abs(a[1] - b[1]) - 1) < tol:
            if abs(a[2] - b[2]) < tol or abs(abs(a[2] - b[2]) - 1) < tol:
                return True

    return False


def apply_symmetry_operations(atoms, spg, tol=1e-4):
    """Figure out the effects of all the symmetry operations."""
    trans = spg._dataset["translations"]
    rots = spg._dataset["rotations"]
    # A = atoms.lattice_mat
    coords = atoms.frac_coords % 1
    nat = coords.shape[0]

    found_everything = True

    permutations = []
    for rot, tran in zip(rots, trans):
        # print("rot", rot)
        # print("tran", tran)
        order = np.zeros(nat, dtype=int)
        found = False
        for at in range(nat):
            pos_new = np.dot(coords[at, :], rot.transpose()) + tran
            pos_new = pos_new % 1

            for at2 in range(nat):
                if check_match(pos_new, coords[at2, :]):
                    found = True
                    order[at] = at2
        if not found:
            found_everything = False
        permutations.append(order)

    return permutations, found_everything


def get_unique_magnetic_structures(
    atoms, supercell_dim=[1, 1, 1], magnetic_ions=None, noferri=True
):
    """
    Generate supercells with unique magnetic orderings.

    noferri=False to get ferrimagnetic configurations.
    """
    if magnetic_ions is None:
        magnetic_ions = set(atoms.elements)

    ss = atoms.make_supercell(dim=supercell_dim)
    spg = Spacegroup3D(ss)

    # apply symmetry with various tolerances until we find one that works
    worked = False
    for tol in [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]:
        permutations, worked = apply_symmetry_operations(ss, spg, tol=tol)
        if worked:
            print("applied sym ", tol)
            break
    if not worked:
        print("error in apply_symmetry_operations")

    print("number of sym permutations:", len(permutations))
    nat = ss.frac_coords.shape[0]

    magnetic_list = []
    magnetic_count = 0
    mag_dict = {}
    for i, el in enumerate(ss.elements):
        if el in magnetic_ions:
            magnetic_list.append(True)
            mag_dict[magnetic_count] = i
            magnetic_count += 1

        else:
            magnetic_list.append(False)

    print("magnetic count: ", magnetic_count)
    if magnetic_count == 0:
        print("no magnetic ions, what are you doing????")
        return

    # generate all magnetic configurations
    magnetic_list = []
    for i in range(2 ** (magnetic_count)):
        magnetic_list.append(np.zeros(magnetic_count))
        # magnetic_list[-1][0] = 5.0

    tmp = "00000000000000000000000000000000000000000000000000000000000000"
    if magnetic_count > 0:
        for i in range(2 ** (magnetic_count)):
            binary_int = bin(i).replace("0b", "")  # convert to binary
            total_int = (
                tmp
                + binary_int
            )

            for ii, d in enumerate(total_int[-magnetic_count:]):
                if d == "0":
                    # print(i, ii, d)
                    magnetic_list[i][ii] = 5.0
                else:
                    # print(i, ii, d)
                    magnetic_list[i][ii] = -5.0

    if (
        noferri
    ):  # get rid if ferrimagnetic configurations, only exact AFM and FM
        newlist = []
        for i in range(2 ** (magnetic_count)):
            if (
                np.abs(np.abs(np.sum(magnetic_list[i])) / 5.0 - magnetic_count)
                < 1e-5
                or np.abs(np.sum(magnetic_list[i])) < 1e-5
            ):
                newlist.append(magnetic_list[i])
        magnetic_list = newlist

    # convert to all atoms in cell
    mag_all = []
    for mag in magnetic_list:
        z = np.zeros(nat)
        for i, m in enumerate(mag):
            z[mag_dict[i]] = m
        mag_all.append(z)

    print("generated, now apply sym opps to find unique")
    # apply symmetry ops
    symm_list = []
    for mag in mag_all:
        already_in_list = False
        for p in permutations:
            mag_new = mag[p]
            for s in symm_list:
                if (
                    np.sum(np.abs(s - mag_new)) < 1e-5
                    or np.sum(np.abs(s + mag_new)) < 1e-5
                ):  # either we found the same config, or same config * -1
                    already_in_list = True
                    break
        if not already_in_list:  # then we need it
            symm_list.append(mag)

    print("number of unique configs: ", len(symm_list))
    return symm_list, ss


"""
if __name__ == "__main__":
    from jarvis.core.atoms import Atoms
    box = [[2.715, 0, 0], [0, 2.715, 0], [0, 0, 2.715]]
    coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
    elements = ["Mn", "Si"]
    MnSi_cube = Atoms(lattice_mat=box, coords=coords, elements=elements)
    from jarvis.db.figshare import get_jid_data

    from jarvis.core.atoms import get_supercell_dims
    atoms = Atoms.from_dict(
        get_jid_data(jid="JVASP-78681", dataset="dft_3d")["atoms"]
    )
    dim = get_supercell_dims(atoms)
    print("dim=", dim)
    dim = [2, 2, 2]
    print("dim=", dim)
    symm_list, ss = get_unique_magnetic_structures(
        atoms, supercell_dim=dim, magnetic_ions=["Mn"]
    )
    print("symm_list", symm_list)
    print("supercell", ss)
"""
