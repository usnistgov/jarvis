from __future__ import division, unicode_literals

"""
Surface abd defect structure generator
See http://iopscience.iop.org/article/10.1088/1361-648X/aadaff/meta
and
See DOI: 10.1016/j.cpc.2015.03.015
"""


import argparse
import os
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.defects.generators import VacancyGenerator
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
from pymatgen.io.vasp import Kpoints

try:
    from pymatgen.ext.matproj import MPRester
except:
    pass


def vac_antisite_def_struct_gen(c_size=15, mpid="", struct=None, write_file=True):
    """
    Vacancy, antisite generator

    Args:
         c_size: cell size
         struct: Structure object or
         mpid: materials project id
    Returns:
            def_str: defect structures in Poscar object format
    """
    def_str = []
    if struct == None:
        with MPRester() as mp:
            struct = mp.get_structure_by_material_id(mpid)
        if mpid == "":
            print("Provide structure")
    c_size = c_size
    prim_struct_sites = len(struct.sites)
    struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    dim1 = int((float(c_size) / float(max(abs(struct.lattice.matrix[0]))))) + 1
    dim2 = int(float(c_size) / float(max(abs(struct.lattice.matrix[1])))) + 1
    dim3 = int(float(c_size) / float(max(abs(struct.lattice.matrix[2])))) + 1
    cellmax = max(dim1, dim2, dim3)
    conv_struct_sites = len(struct.sites)
    conv_prim_rat = int(conv_struct_sites / prim_struct_sites)
    sc_scale = [dim1, dim2, dim3]
    print("sc_scale", sc_scale)

    tmp = struct.copy()
    tmp.make_supercell(sc_scale)
    sc_tmp = tmp  # Poscar(tmp).structure .make_supercell(list(sc_scale))
    scs = list(VacancyGenerator(struct))
    supercell = Poscar(sc_tmp)
    supercell.comment = str("bulk") + str("@") + str("cellmax") + str(cellmax)
    def_str.append(supercell)
    if write_file == True:
        supercell.write_file("POSCAR-" + str("bulk") + str(".vasp"))

    for i in range(len(scs)):
        sc = scs[i].generate_defect_structure(sc_scale)
        poscar = Poscar(sc)  # mpvis.get_poscar(sc)
        pmg_name = str(scs[i].name).split("_")
        sitespecie = pmg_name[1]
        mult = pmg_name[2].split("mult")[1]
        name = (
            str("vacancy_")
            + str(i + 1)
            + str("_mult-")
            + str(mult)
            + str("_sitespecie-")
            + str(sitespecie)
            + str("@cellmax")
            + str(cellmax)
        )
        poscar.comment = str(name)
        def_str.append(poscar)
        if write_file == True:
            filename = (
                str("POSCAR-")
                + str("vacancy_")
                + str(i + 1)
                + str("_mult-")
                + str(mult)
                + str("_sitespecie-")
                + str(sitespecie)
                + str(".vasp")
            )
            poscar.write_file(filename)

    return def_str


def pmg_surfer(
    mpid="", vacuum=15, mat=None, max_index=1, min_slab_size=15, write_file=True
):
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
        with MPRester() as mp:
            mat = mp.get_structure_by_material_id(mpid)
        if mpid == "":
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


def surfer(mpid="", vacuum=15, layers=2, mat=None, max_index=1, write_file=True):
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
        with MPRester() as mp:
            mat = mp.get_structure_by_material_id(mpid)
        if mpid == "":
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


def main():
    pp = vac_antisite_def_struct_gen(cellmax=2, mpid="mp-134")
    # pp=vac_intl(cellmax=128,mpid='mp-134')


#    ss=surfer(mpid='mp-134')
#    print ss
# main()
# strt = Structure.from_file("POSCAR")
# vac_antisite_def_struct_gen(struct=strt)
