from __future__ import division, unicode_literals

"""
A script with tools for computing point defect concentrations.
Manual and citation for the script, DOI: 10.1016/j.cpc.2015.03.015
"""


import argparse
import os
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp import Poscar
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.defects.core import Vacancy,Interstitial
from pymatgen.analysis.local_env import ValenceIonicRadiusEvaluator
#from pymatgen.analysis.defects.point_defects import Vacancy,Interstitial,ValenceIonicRadiusEvaluator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import glob
from pymatgen.io.vasp import Poscar
from monty.serialization import loadfn, dumpfn
from monty.json import MontyEncoder, MontyDecoder
#from pymatgen.analysis.defects.point_defects import Vacancy
#from pymatgen.io.vasp.sets import MPGGAVaspInputSet
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



def get_sc_scale(inp_struct, final_site_no):
    lengths = inp_struct.lattice.abc
    no_sites = inp_struct.num_sites
    mult = (final_site_no/no_sites*lengths[0]*lengths[1]*lengths[2]) ** (1/3)
    num_mult = [int(round(mult/l)) for l in lengths]
    num_mult = [i if i > 0 else 1 for i in num_mult]
    return num_mult


def vac_antisite_def_struct_gen(c_size=15,mpid='',struct=None):
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
    #print ("in vac_def cellmax=",cell
    prim_struct_sites = len(struct.sites)
    struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    conv_struct_sites = len(struct.sites)
    conv_prim_rat = int(conv_struct_sites/prim_struct_sites)
    #sc_scale = get_sc_scale(struct,cellmax)
    sc_scale=[dim1,dim2,dim3]
    #sc_scale=[cellmax,cellmax,cellmax]
    print ("sc_scale",sc_scale)
    #mpvis = MPGGAVaspInputSet()

    # Begin defaults: All default settings.
    #blk_vasp_incar_param = {'IBRION':-1,'EDIFF':1e-4,'EDIFFG':0.001,'NSW':0,}
    #def_vasp_incar_param = {'ISIF':2,'NELM':99,'IBRION':2,'EDIFF':1e-6, 
    #                        'EDIFFG':0.001,'NSW':40,}
    #kpoint_den = 6000
    # End defaults
    
    #ptcr_flag = True
    #try:
    #    potcar = mpvis.get_potcar(struct)
    #except:
    #    print ("VASP POTCAR folder not detected.\n" \
    #          "Only INCAR, POSCAR, KPOINTS are generated.\n" \
    #          "If you have VASP installed on this system, \n" \
    #          "refer to pymatgen documentation for configuring the settings.")
    #    ptcr_flag = False





    struct_valrad_eval = ValenceIonicRadiusEvaluator(struct)
    val = struct_valrad_eval.valences
    rad = struct_valrad_eval.radii
    struct_val = val
    struct_rad = rad




    vac = Vacancy(struct, {}, {})
    scs = vac.make_supercells_with_defects(sc_scale)
    #site_no = scs[0].num_sites
    #if site_no > cellmax:
    #    max_sc_dim = max(sc_scale)
    #    i = sc_scale.index(max_sc_dim)
    #    sc_scale[i] -= 1
    #    scs = vac.make_supercells_with_defects(sc_scale)


    for i in range(len(scs)):
        sc = scs[i]
        poscar = mpvis.get_poscar(sc)
        kpoints = Kpoints.automatic_density(sc,kpoint_den)
        incar = mpvis.get_incar(sc)
        #if ptcr_flag:
        #    potcar = mpvis.get_potcar(sc)

        interdir = mpid
        if not i:
            fin_dir = os.path.join(interdir,'bulk')
            #try:
            #    os.makedirs(fin_dir)
            #except:
            #    pass
            #incar.update(blk_vasp_incar_param)
            #incar.write_file(os.path.join(fin_dir,'INCAR'))
            #poscar.write_file(os.path.join(fin_dir,'POSCAR'))
            poscar.comment=str('bulk')+str('@')+str('cellmax')+str(cellmax)
            def_str.append(poscar)
            poscar.write_file('POSCAR-'+str('bulk')+str(".vasp"))
            #if ptcr_flag:
            #    potcar.write_file(os.path.join(fin_dir,'POTCAR'))
            #kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))
        else:
            blk_str_sites = set(scs[0].sites)
            vac_str_sites = set(sc.sites)
            vac_sites = blk_str_sites - vac_str_sites
            vac_site = list(vac_sites)[0]
            site_mult = int(vac.get_defectsite_multiplicity(i-1)/conv_prim_rat)
            #try:
            #   site_mult = int(vac.get_defectsite_multiplicity(i-1)/conv_prim_rat)
            #except:
            #   site_mult=1
            #   pass
            vac_site_specie = vac_site.specie
            vac_symbol = vac_site.specie.symbol

            vac_dir ='vacancy_{}_mult-{}_sitespecie-{}'.format(str(i),
                    site_mult, vac_symbol)
            fin_dir = os.path.join(interdir,vac_dir)
            #try:
            #    os.makedirs(fin_dir)
            #except:
            #    pass
            #incar.update(def_vasp_incar_param)
            try:
                poscar.comment=str(vac_dir)+str('@')+str('cellmax')+str(cellmax)
            except:
                pass
            pos=poscar
            #pos=poscar.structure
            def_str.append(pos)
            #poscar.write_file(os.path.join(fin_dir,'POSCAR'))
            poscar.write_file('POSCAR-'+str(vac_dir)+str(".vasp"))
            #incar.write_file(os.path.join(fin_dir,'INCAR'))
            #if ptcr_flag:
            #    potcar.write_file(os.path.join(fin_dir,'POTCAR'))
            #kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))

            # Antisite generation at all vacancy sites
            struct_species = scs[0].types_of_specie
            for specie in set(struct_species)-set([vac_site_specie]):
                subspecie_symbol = specie.symbol
                anti_struct = sc.copy()
                anti_struct.append(specie, vac_site.frac_coords)
                poscar = mpvis.get_poscar(anti_struct)
                incar = mpvis.get_incar(anti_struct)
                incar.update(def_vasp_incar_param)
                as_dir ='antisite_{}_mult-{}_sitespecie-{}_subspecie-{}'.format(
                        str(i), site_mult, vac_symbol, subspecie_symbol)
                fin_dir = os.path.join(interdir,as_dir)
                #try:
                #    os.makedirs(fin_dir)
                #except:
                #    pass
                poscar.comment=str(as_dir)+str('@')+str('cellmax')+str(cellmax)
                pos=poscar
                #pos=poscar.structure
                def_str.append(pos)
                #poscar.write_file(os.path.join(fin_dir,'POSCAR'))
                poscar.write_file('POSCAR-'+str(as_dir)+str(".vasp"))
                #incar.write_file(os.path.join(fin_dir,'INCAR'))
                #if ptcr_flag:
                #        potcar.write_file(os.path.join(fin_dir,'POTCAR'))
                #kpoints.write_file(os.path.join(fin_dir,'KPOINTS'))
    #try:
    #    struct.make_supercell(sc_scale) 
    #    intl = Interstitial(struct, val, rad)
    #    cell_arr=sc_scale
    #    for el in struct.composition.elements:
    #        scs = intl.make_supercells_with_defects(1,el)
    #        #scs = intl.make_supercells_with_defects(cell_arr,el)
    #        for i in range(1,len(scs)):
    #           pos=Poscar(scs[i])
    #           pos.comment=str('intl_')+str('cellmax')+str(cellmax)+str('@')+str(intl.get_defectsite_coordination_number(i-1))+str('Element')+str(el)
    #           def_str.append(pos)
    #           pos.write_file('POSCAR-'+str('intl_')+str('cellmax')+str(cellmax)+str('@')+str(intl.get_defectsite_coordination_number(i-1))+str('Element')+str(el))

    #except:
    #      pass

    return def_str
def pmg_surfer(mpid='',vacuum=15,mat=None,max_index=1,min_slab_size=15):
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
        pos.write_file(filename=str('POSCAR-')+str("Surf-")+str(surf_name)+str('.vasp'))
        structures.append(pos)

    return structures
def surfer(mpid='',vacuum=15,layers=2,mat=None,max_index=1):
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
        pos.write_file(filename=str('POSCAR-')+str("Surf-")+str(surf_name)+str('.vasp'))
        structures.append(pos)

    return structures


def vac_intl(cellmax=2,mpid='',struct=None):


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

def main():
    pp=vac_antisite_def_struct_gen(cellmax=2,mpid='mp-134')
    #pp=vac_intl(cellmax=128,mpid='mp-134')
#    ss=surfer(mpid='mp-134')
#    print ss
#main()
#strt = Structure.from_file("POSCAR")
#vac_antisite_def_struct_gen(struct=strt)

