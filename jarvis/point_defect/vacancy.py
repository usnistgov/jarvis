from __future__ import division, unicode_literals

"""
Surface and defect structure generator
See http://iopscience.iop.org/article/10.1088/1361-648X/aadaff/meta
and
See DOI: 10.1016/j.cpc.2015.03.015
"""


import argparse
import os
from pymatgen.analysis.defects.generators import VacancyGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import glob
from pymatgen.io.vasp import Poscar
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from monty.serialization import loadfn, MontyDecoder


def chempot_struct(element="Al"):
    chempot_json_file = str(os.path.join(os.path.dirname(__file__), "chem_pot.json"))
    chempot_json = loadfn(chempot_json_file, cls=MontyDecoder)
    id, struct = chempot_json[element]
    return id, struct


def vac_antisite_def_struct_gen(c_size=15, struct=None, write_file=True):
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


if __name__ == "__main__":
    pos = str(
        os.path.join(
            os.path.dirname(__file__),
            "../vasp/examples/SiOptb88/MAIN-MBJ-bulk@mp_149/POSCAR",
        )
    )
    mat = Structure.from_file(pos)
    x = vac_antisite_def_struct_gen(struct=mat, write_file=False)
    print(x)
    id, s = chempot_struct("Al")
    print(id, s)
