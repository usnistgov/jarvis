"""Module for setting MAGMOM and AFM/FM orderings."""

from jarvis.analysis.structure.spacegroup import Spacegroup3D
import numpy as np
from collections import defaultdict
from jarvis.core.utils import check_match


class MagneticOrdering(object):
    """Provide modules for enumerating magnetic ordering analysis."""

    def __init__(self, atoms=None):
        """Initialize with Atoms."""
        self.atoms = atoms

    def apply_symmetry_operations(self, atoms, spg, tol=1e-4):
        """Figure out the effects of all the symmetry operations."""
        trans = spg._dataset["translations"]
        rots = spg._dataset["rotations"]
        coords = atoms.frac_coords % 1
        nat = atoms.num_atoms
        found_everything = True
        permutations = []
        for rot, tran in zip(rots, trans):
            order = np.zeros(nat, dtype=int)
            found = False
            for at in range(nat):
                pos_new = np.dot(coords[at, :], rot.transpose()) + tran
                pos_new = pos_new % 1
                for at2 in range(nat):
                    if check_match(pos_new, coords[at2, :], tol=tol):
                        found = True
                        order[at] = at2
            if not found:
                found_everything = False
            permutations.append(order)
        return permutations, found_everything

    def get_unique_magnetic_structures(
        self,
        atoms,
        supercell_dim=[1, 1, 1],
        magnetic_ions=None,
        noferri=True,
        magmom=3.0,
    ):
        """
        Generate supercells with unique magnetic orderings.

        noferri=False to get ferrimagnetic configurations.
        """
        if magnetic_ions is None:
            magnetic_ions = set(atoms.elements)

        ss = atoms.make_supercell(dim=supercell_dim)
        # spg = Spacegroup3D(atoms)
        spg = Spacegroup3D(atoms)  # kfg

        # Apply symmetry with various tolerances until we find one that works
        worked = False
        for tol in [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]:
            permutations, worked = self.apply_symmetry_operations(
                ss, spg, tol=tol
            )
            if worked:
                print("applied sym ", tol)
                break
        if not worked:
            print("error in apply_symmetry_operations")

        print("number of sym permutations:", len(permutations))
        nat = ss.num_atoms

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
                total_int = tmp + binary_int

                for ii, d in enumerate(total_int[-magnetic_count:]):
                    if d == "0":
                        # print(i, ii, d)
                        magnetic_list[i][ii] = magmom
                    else:
                        # print(i, ii, d)
                        magnetic_list[i][ii] = -1 * magmom

        if (
            noferri
        ):  # get rid if ferrimagnetic configurations, only exact AFM and FM
            newlist = []
            for i in range(2 ** (magnetic_count)):
                if (
                    np.abs(
                        np.abs(np.sum(magnetic_list[i])) / abs(magmom)
                        - magnetic_count
                    )
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

    def tc_mean_field(self, energies=[-2, -1]):
        """Curie temperature using mean-field theory."""
        info = {}
        mag_atoms = self.get_mag_ions()
        deltaE = max(energies) - min(energies)
        print(deltaE)
        kB = 8.617333262e-05
        elements_dict = defaultdict(int)
        for i in self.atoms.elements:
            if i in mag_atoms:
                elements_dict[i] += 1
        n_mag_elements = sum(elements_dict.values())
        Tc = deltaE / (3 * kB) / n_mag_elements
        # Tc = 2 * deltaE / (3 * kB) / n_mag_elements
        info["Tc"] = Tc
        info["deltaE"] = deltaE
        info["n_mag_elements"] = n_mag_elements
        return info

    def get_mag_ions(self):
        """List all magnetic atoms in the Atoms object."""
        all_mag_elements = [
            "Ti",
            "V",
            "Cr",
            "Mn",
            "Fe",
            "Co",
            "Ni",
            "Cu",
            "Ru",
            "Ir",
            "Rh",
            "Os",
            "Rb",
            "Sc",
        ]
        els = self.atoms.elements
        mag_ions = list(set(all_mag_elements).intersection(set(els)))
        return mag_ions

    def get_minimum_configs(self, min_configs=3, enforce_primitive=True):
        """Get minimum number of spin structures for Tc calculations."""
        atoms = self.atoms
        if enforce_primitive:
            atoms = atoms.get_primitive_atoms
        lengths = atoms.lattice.lat_lengths()
        index_to_expand = np.argsort(lengths)
        dim = np.array([1, 1, 1])
        mag_ions = self.get_mag_ions()
        symm_list, ss = self.get_unique_magnetic_structures(
            atoms, supercell_dim=dim, magnetic_ions=mag_ions, noferri=True
        )
        # kfg
        if len(symm_list) < min_configs:
            symm_list, ss = self.get_unique_magnetic_structures(
                atoms, supercell_dim=dim, magnetic_ions=mag_ions, noferri=True
            )

        count = 0
        while len(symm_list) < min_configs:
            dim[index_to_expand[count]] += 1
            symm_list, ss = self.get_unique_magnetic_structures(
                atoms,
                supercell_dim=dim,
                magnetic_ions=self.get_mag_ions(),
                noferri=True,
            )
            # kfg
            if len(symm_list) < min_configs:
                symm_list, ss = self.get_unique_magnetic_structures(
                    atoms,
                    supercell_dim=dim,
                    magnetic_ions=self.get_mag_ions(),
                    noferri=False,
                )

            count = count + 1
            if count > 2:
                count = 0
            print("Supercell dimension", dim)
        return symm_list, ss
