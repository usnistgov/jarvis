"""Module to prepare input file for prismatic."""

# Use jarvis.core.atoms.crop_squre for non rectangles


def write_prismatic_xyz(atoms=None):
    """Write Atoms in prismatic format."""
    from jarvis.core.specie import Specie

    line = str(atoms.num_atoms) + "\n"
    abc = atoms.lattice.abc
    line += str(abc[0]) + " " + str(abc[1]) + " " + str(abc[2]) + "\n"
    for i, j in zip(atoms.elements, atoms.cart_coords):
        line += (
            str(Specie(i).Z)
            + " "
            + str(j[0])
            + " "
            + str(j[1])
            + " "
            + str(j[2])
            + str("  1 .076\n")
        )
    line += "-1\n"
    return line
