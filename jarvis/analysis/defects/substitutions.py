"""Module to substitute atoms."""

from jarvis.core.atoms import Atoms
def substitute_atoms(atoms=None, element="Na", site=0):
    """Substitute element in Atoms class."""
    elements = atoms.elements
    elements[site] = element
    return Atoms(
        lattice_mat=atoms.lattice_mat,
        elements=elements,
        coords=atoms.coords,
        cartesian=atoms.cartesian,
    )
