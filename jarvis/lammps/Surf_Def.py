from __future__ import division, unicode_literals

"""
Surface abd defect structure generator
See DOI: 10.1016/j.cpc.2015.03.015
"""


import argparse
import os
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.defects.point_defects import ValenceIonicRadiusEvaluator,Vacancy,Interstitial
#from pymatgen.analysis.local_env import ValenceIonicRadiusEvaluator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import glob
from pymatgen.io.vasp import Poscar
from monty.serialization import loadfn, dumpfn
from monty.json import MontyEncoder, MontyDecoder
#from pymatgen.analysis.defects.point_defects import Vacancy
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.lattice.surface import surface
from pymatgen.core.surface import  Slab, SlabGenerator, generate_all_slabs,get_symmetrically_distinct_miller_indices
from pymatgen.io.vasp import Kpoints
from pymatgen.io.vasp import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.defects.dilute_solution_model import \
            compute_defect_density, solute_defect_density
try:
  from pymatgen.analysis.defects.point_defects import Vacancy,Interstitial,ValenceIonicRadiusEvaluator
except:
  from pymatgen.analysis.defects.core import Vacancy,Interstitial
  pass

try: 
  from pymatgen.ext.matproj import MPRester
except:
    pass

def get_sc_scale(inp_struct, final_site_no):
    """
    scale up structure

    Args:
        inp_struct: initial input structure
        final_site_no: final number of atoms
    Returns:
           number of multipliers
    """
    lengths = inp_struct.lattice.abc
    no_sites = inp_struct.num_sites
    mult = (final_site_no/no_sites*lengths[0]*lengths[1]*lengths[2]) ** (1/3)
    num_mult = [int(round(mult/l)) for l in lengths]
    num_mult = [i if i > 0 else 1 for i in num_mult]
    return num_mult


def vac_antisite_def_struct_gen(c_size=15,mpid='',struct=None,write_file=True):
    """
    Vacancy, antisite generator

    Args:
         c_size: cell size
         struct: Structure object or
         mpid: materials project id
    Returns:
            def_str: defect structures in Poscar object format
    """     
    def_str=[]
    if struct ==None:
        with MPRester() as mp:
                struct = mp.get_structure_by_material_id(mpid)
        if mpid == '':
           print ("Provide structure")
    c_size=c_size
    dim1=int((float(c_size)/float( max(abs(struct.lattice.matrix[0])))))+1
    dim2=int(float(c_size)/float( max(abs(struct.lattice.matrix[1]))))+1
    dim3=int(float(c_size)/float( max(abs(struct.lattice.matrix[2]))))+1
    cellmax=max(dim1,dim2,dim3)
    prim_struct_sites = len(struct.sites)
    struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    conv_struct_sites = len(struct.sites)
    conv_prim_rat = int(conv_struct_sites/prim_struct_sites)
    sc_scale=[dim1,dim2,dim3]
    print ("sc_scale",sc_scale)

    struct_valrad_eval = ValenceIonicRadiusEvaluator(struct)
    val = struct_valrad_eval.valences
    rad = struct_valrad_eval.radii
    struct_val = val
    struct_rad = rad

    vac = Vacancy(struct, {}, {})
    scs = vac.make_supercells_with_defects(sc_scale)


    for i in range(len(scs)):
        sc = scs[i]
        poscar = Poscar(sc) #mpvis.get_poscar(sc)

        interdir = mpid
        if not i:
            fin_dir = os.path.join(interdir,'bulk')
            poscar.comment=str('bulk')+str('@')+str('cellmax')+str(cellmax)
            def_str.append(poscar)
            if write_file==True:
                poscar.write_file('POSCAR-'+str('bulk')+str(".vasp"))
        else:
            blk_str_sites = set(scs[0].sites)
            vac_str_sites = set(sc.sites)
            vac_sites = blk_str_sites - vac_str_sites
            vac_site = list(vac_sites)[0]
            site_mult = int(vac.get_defectsite_multiplicity(i-1)/conv_prim_rat)
            vac_site_specie = vac_site.specie
            vac_symbol = vac_site.specie.symbol

            vac_dir ='vacancy_{}_mult-{}_sitespecie-{}'.format(str(i),
                    site_mult, vac_symbol)
            fin_dir = os.path.join(interdir,vac_dir)
            try:
                poscar.comment=str(vac_dir)+str('@')+str('cellmax')+str(cellmax)
            except:
                pass
            pos=poscar
            def_str.append(pos)
            if write_file==True:
                poscar.write_file('POSCAR-'+str(vac_dir)+str(".vasp"))
            struct_species = scs[0].types_of_specie
            for specie in set(struct_species)-set([vac_site_specie]):
                subspecie_symbol = specie.symbol
                anti_struct = sc.copy()
                anti_struct.append(specie, vac_site.frac_coords)
                
                poscar = Poscar(anti_struct)
                as_dir ='antisite_{}_mult-{}_sitespecie-{}_subspecie-{}'.format(
                        str(i), site_mult, vac_symbol, subspecie_symbol)
                fin_dir = os.path.join(interdir,as_dir)
                poscar.comment=str(as_dir)+str('@')+str('cellmax')+str(cellmax)
                pos=poscar
                def_str.append(pos)
                if write_file==True:
                    poscar.write_file('POSCAR-'+str(as_dir)+str(".vasp"))

    return def_str

def pmg_surfer(mpid='',vacuum=15,mat=None,max_index=1,min_slab_size=15,write_file=True):
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
        if mpid == '':
          print ('Provide structure')

    sg_mat = SpacegroupAnalyzer(mat)
    mat_cvn = sg_mat.get_conventional_standard_structure()
    mat_cvn.sort()
    indices = get_symmetrically_distinct_miller_indices(mat_cvn, max_index)
    #ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)

    structures=[]
    pos=Poscar(mat_cvn)
    try:
       pos.comment=str('sbulk')+str('@')+str('vac')+str(vacuum)+str('@')+str('size')+str(min_slab_size)
    except:
       pass
    structures.append(pos)
    if write_file==True:
        mat_cvn.to(fmt='poscar',filename=str('POSCAR-')+str('cvn')+str('.vasp'))
    for i in indices:
        slab=SlabGenerator(initial_structure = mat_cvn, miller_index=i, min_slab_size= min_slab_size, min_vacuum_size=vacuum , lll_reduce=False, center_slab=True, primitive=False).get_slab()
        normal_slab = slab.get_orthogonal_c_slab()
        slab_pymatgen = Poscar(normal_slab).structure
        #ase_slab.center(vacuum=vacuum, axis=2)
        #slab_pymatgen = AseAtomsAdaptor().get_structure(ase_slab)
        xy_size=min_slab_size
        dim1=int((float(xy_size)/float( max(abs(slab_pymatgen.lattice.matrix[0])))))+1
        dim2=int(float(xy_size)/float( max(abs(slab_pymatgen.lattice.matrix[1]))))+1
        slab_pymatgen.make_supercell([dim1,dim2,1])
        slab_pymatgen.sort()
        surf_name='_'.join(map(str,i))
        pos=Poscar(slab_pymatgen)
        try:
           pos.comment=str("Surf-")+str(surf_name)+str('@')+str('vac')+str(vacuum)+str('@')+str('size')+str(min_slab_size)
        except:
           pass
        if write_file==True:
             pos.write_file(filename=str('POSCAR-')+str("Surf-")+str(surf_name)+str('.vasp'))
        structures.append(pos)

    return structures

def surfer(mpid='',vacuum=15,layers=2,mat=None,max_index=1,write_file=True):
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
        if mpid == '':
          print ('Provide structure')

    sg_mat = SpacegroupAnalyzer(mat)
    mat_cvn = sg_mat.get_conventional_standard_structure()
    mat_cvn.sort()
    indices = get_symmetrically_distinct_miller_indices(mat_cvn, max_index)
    ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)

    structures=[]
    pos=Poscar(mat_cvn)
    try:
       pos.comment=str('sbulk')+str('@')+str('vac')+str(vacuum)+str('@')+str('layers')+str(layers)
    except:
       pass
    structures.append(pos)
    if write_file==True:
        mat_cvn.to(fmt='poscar',filename=str('POSCAR-')+str('cvn')+str('.vasp'))
    for i in indices:
        ase_slab = surface(ase_atoms, i, layers)
        ase_slab.center(vacuum=vacuum, axis=2)
        slab_pymatgen = AseAtomsAdaptor().get_structure(ase_slab)
        slab_pymatgen.sort()
        surf_name='_'.join(map(str,i))
        pos=Poscar(slab_pymatgen)
        try:
           pos.comment=str("Surf-")+str(surf_name)+str('@')+str('vac')+str(vacuum)+str('@')+str('layers')+str(layers)
        except:
           pass
        if write_file==True:
            pos.write_file(filename=str('POSCAR-')+str("Surf-")+str(surf_name)+str('.vasp'))
        structures.append(pos)

    return structures


def vac_intl(cellmax=2,mpid='',struct=None):

    """
    Vacancy and interstitial generator

    Args:
        cellmax: maximum cell size
        struct: Structure object
    Returns:
            def_str: defect structures
    """

    if struct ==None:
        with MPRester() as mp:
                struct = mp.get_structure_by_material_id(mpid)
        if mpid == '':
           print ("Provide structure")
    sg_mat = SpacegroupAnalyzer(struct)
    struct = sg_mat.get_conventional_standard_structure()
    def_str=[]
    count=0
    cell_arr=[cellmax,cellmax,cellmax]
    vac = Vacancy(struct, {}, {})
    for el in list(struct.symbol_set):
   
        scs = vac.make_supercells_with_defects(cell_arr,el)
        for i in range(len(scs)):
            if i==0:
               pos=Poscar(scs[i])
               pos.comment=str('bulk')+str('.')+str('cellmax')+str(cellmax)
               if count==0:
                   def_str.append(pos)
                   count=count+1
            else:
               pos=Poscar(scs[i])
               pos.comment=str('vac')+str('cellmax')+str(cellmax)+str('@')+str(vac.get_defectsite_multiplicity(i))+str('Element')+str(el)
               if pos not in def_str:
                   def_str.append(pos)
    struct_valrad_eval = ValenceIonicRadiusEvaluator(struct)
    val = struct_valrad_eval.valences
    rad = struct_valrad_eval.radii
    struct_val = val
    struct_rad = rad
    intl = Interstitial(struct, val, rad)

    for el in struct.composition.elements:
        scs = intl.make_supercells_with_defects(cell_arr,el)
        for i in range(1,len(scs)):
           pos=Poscar(scs[i])
           pos.comment=str('intl')+str('cellmax')+str(cellmax)+str('@')+str(intl.get_defectsite_coordination_number(i-1))+str('Element')+str(el)
           if pos not in def_str:
                   def_str.append(pos)
    print (len(def_str))
    return def_str
def main():
    pp=vac_antisite_def_struct_gen(cellmax=2,mpid='mp-134')
    #pp=vac_intl(cellmax=128,mpid='mp-134')
#    ss=surfer(mpid='mp-134')
#    print ss
#main()
#strt = Structure.from_file("POSCAR")
#vac_antisite_def_struct_gen(struct=strt)

