"""Modules for making point-defect substituions."""
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.core.utils import rand_select
from jarvis.core.atoms import Atoms, get_supercell_dims


def generate_defect(
    atoms=None,
    enforce_c_size=10.0,
    on_conventional_cell=False,
    extend=1,
    subs_element="Al",
    selected_element=None,
):
    """Provide function to generate substitution defects."""
    if on_conventional_cell:
        atoms = Spacegroup3D(atoms).conventional_standard_structure
    if enforce_c_size is not None:
        dims = get_supercell_dims(
            atoms, enforce_c_size=enforce_c_size, extend=extend
        )
        supercell_size = [dims[0], dims[1], dims[2]]
        # atoms = atoms.make_supercell(supercell_size)
    spg = Spacegroup3D(atoms)
    wyckoffs = spg._dataset["wyckoffs"]
    atoms.props = wyckoffs
    supercell = atoms.make_supercell(supercell_size)
    props = rand_select(supercell.props)
    subs = []
    # print(props)
    for i, j in props.items():
        info = {}
        elements = supercell.elements.copy()
        if selected_element is not None:
            if elements[j] == selected_element:

                elements[j] = subs_element

                a = Atoms(
                    lattice_mat=supercell.lattice_mat,
                    coords=supercell.coords,
                    cartesian=supercell.cartesian,
                    elements=elements,
                )
                info["props"] = props
                info["atoms"] = atoms.to_dict()
                info["defect_atoms"] = a.to_dict()
                info["supercell_size"] = list(supercell_size)
                # print(a.elements)
                subs.append(info)
        else:
            elements[j] = subs_element
            a = Atoms(
                lattice_mat=supercell.lattice_mat,
                coords=supercell.coords,
                cartesian=supercell.cartesian,
                elements=elements,
            )
            info["props"] = props
            info["atoms"] = atoms.to_dict()
            info["defect_atoms"] = a.to_dict()
            info["supercell_size"] = list(supercell_size)
            # print(a.elements)
            subs.append(info)
    return subs


"""
x = generate_defect(atoms=a, selected_element="Br")
"""
