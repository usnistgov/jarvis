from __future__ import division, unicode_literals

"""
Surface and defect structure generator
See http://iopscience.iop.org/article/10.1088/1361-648X/aadaff/meta
and
See DOI: 10.1016/j.cpc.2015.03.015
"""


import argparse
import os
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import glob
from pymatgen.io.vasp import Poscar
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.lattice.surface import surface
from pymatgen.core.surface import (
    Slab,
    SlabGenerator,
    generate_all_slabs,
    get_symmetrically_distinct_miller_indices,
)


def pmg_surfer(vacuum=15, mat=None, max_index=1, min_slab_size=15, write_file=True):
    """
    Pymatgen surface builder for a Poscar

    Args:
        vacuum: vacuum region
        mat: Structure object
        max_index: maximum miller index
        min_slab_size: minimum slab size

    Returns:
           structures: list of surface Structure objects
    """

    if mat == None:
        print("Provide structure")

    sg_mat = SpacegroupAnalyzer(mat)
    mat_cvn = sg_mat.get_conventional_standard_structure()
    mat_cvn.sort()
    indices = get_symmetrically_distinct_miller_indices(mat_cvn, max_index)

    structures = []
    pos = Poscar(mat_cvn)
    try:
        pos.comment = (
            str("sbulk")
            + str("@")
            + str("vac")
            + str(vacuum)
            + str("@")
            + str("size")
            + str(min_slab_size)
        )
    except:
        pass
    structures.append(pos)
    if write_file == True:
        mat_cvn.to(fmt="poscar", filename=str("POSCAR-") + str("cvn") + str(".vasp"))
    for i in indices:
        slab = SlabGenerator(
            initial_structure=mat_cvn,
            miller_index=i,
            min_slab_size=min_slab_size,
            min_vacuum_size=vacuum,
            lll_reduce=False,
            center_slab=True,
            primitive=False,
        ).get_slab()
        normal_slab = slab.get_orthogonal_c_slab()
        slab_pymatgen = Poscar(normal_slab).structure
        xy_size = min_slab_size
        dim1 = (
            int((float(xy_size) / float(max(abs(slab_pymatgen.lattice.matrix[0]))))) + 1
        )
        dim2 = (
            int(float(xy_size) / float(max(abs(slab_pymatgen.lattice.matrix[1])))) + 1
        )
        slab_pymatgen.make_supercell([dim1, dim2, 1])
        slab_pymatgen.sort()
        surf_name = "_".join(map(str, i))
        pos = Poscar(slab_pymatgen)
        try:
            pos.comment = (
                str("Surf-")
                + str(surf_name)
                + str("@")
                + str("vac")
                + str(vacuum)
                + str("@")
                + str("size")
                + str(min_slab_size)
            )
        except:
            pass
        if write_file == True:
            pos.write_file(
                filename=str("POSCAR-") + str("Surf-") + str(surf_name) + str(".vasp")
            )
        structures.append(pos)

    return structures


def surfer(vacuum=15, layers=2, mat=None, max_index=1, write_file=True):
    """
    ASE surface bulder

    Args:
        vacuum: vacuum region
        mat: Structure object
        max_index: maximum miller index
        min_slab_size: minimum slab size

    Returns:
           structures: list of surface Structure objects
    """

    if mat == None:
        print("Provide structure")

    sg_mat = SpacegroupAnalyzer(mat)
    mat_cvn = sg_mat.get_conventional_standard_structure()
    mat_cvn.sort()
    indices = get_symmetrically_distinct_miller_indices(mat_cvn, max_index)
    ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)

    structures = []
    pos = Poscar(mat_cvn)
    try:
        pos.comment = (
            str("sbulk")
            + str("@")
            + str("vac")
            + str(vacuum)
            + str("@")
            + str("layers")
            + str(layers)
        )
    except:
        pass
    structures.append(pos)
    if write_file == True:
        mat_cvn.to(fmt="poscar", filename=str("POSCAR-") + str("cvn") + str(".vasp"))
    for i in indices:
        ase_slab = surface(ase_atoms, i, layers)
        ase_slab.center(vacuum=vacuum, axis=2)
        slab_pymatgen = AseAtomsAdaptor().get_structure(ase_slab)
        slab_pymatgen.sort()
        surf_name = "_".join(map(str, i))
        pos = Poscar(slab_pymatgen)
        try:
            pos.comment = (
                str("Surf-")
                + str(surf_name)
                + str("@")
                + str("vac")
                + str(vacuum)
                + str("@")
                + str("layers")
                + str(layers)
            )
        except:
            pass
        if write_file == True:
            pos.write_file(
                filename=str("POSCAR-") + str("Surf-") + str(surf_name) + str(".vasp")
            )
        structures.append(pos)

    return structures


if __name__ == "__main__":
    pos = str(
        os.path.join(
            os.path.dirname(__file__),
            "../vasp/examples/SiOptb88/MAIN-MBJ-bulk@mp_149/POSCAR",
        )
    )
    mat = Structure.from_file(pos)
    x = pmg_surfer(mat=mat, write_file=False)
    print(x)
    y = surfer(mat=mat, write_file=False)
