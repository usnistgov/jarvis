from numpy import matrix
import time
import numpy as np
import glob
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.core.surface import  Slab, SlabGenerator, generate_all_slabs,get_symmetrically_distinct_miller_indices
from ase.lattice.surface import surface
from pymatgen.ext.matproj import MPRester
#from pymatgen.matproj.rest import MPRester
import operator
from pymatgen.core.lattice import Lattice
#from mpinterfaces.MP_lammps import MPINTLammps
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from  jarvis.lammps.jlammps import vac_antisite_def_struct_gen,surfer
import numpy as np,time,json
import sys,os,subprocess,socket
from pymatgen.io.ase import AseAtomsAdaptor
from ase.calculators.lammpsrun import LAMMPS, Prism
import sys,zipfile
import fortranformat as fform
from pymatgen.core.structure import Structure
#from mpinterfaces.MP_lammps import MPINTLammps
import os
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
import sys,shutil
import urllib2
import urllib
from BeautifulSoup import BeautifulSoup
import json
MAPI_KEY = os.environ.get("MAPI_KEY", "")
from pymatgen.matproj.rest import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from math import ceil
from collections import OrderedDict

import matplotlib
matplotlib.use('Agg')

from pymatgen.core.composition import Composition
from pymatgen.phasediagram.entries import PDEntry
#from pymatgen.phasediagram.pdmaker import PhaseDiagram
#from pymatgen.phasediagram.plotter import PDPlotter

#from mpinterfaces import get_struct_from_mp
#from mpinterfaces.MP_lammps import CalibrateLammps
#from mpinterfaces.utils import *
import glob

###################MAIN PART############################################
for file in glob.glob("*.fs"):
    try:
        folder1=str(os.getcwd())+str("/")+str(file)+str("_nist")
        if not os.path.exists(folder1):
           os.makedirs(str(folder1))
        cwd1=str(os.getcwd())
        print ("folder1=",folder1)
        ff=str(file)
        element_ff=[]
        f=open(ff,'r')
        os.chdir(str(folder1))
        list_el=[]
        lines=f.readlines()
        content=(lines[3]).split(" ")
    #content=(lines[3]).split("' '|\n|\r\n")
        for val in content:
        
            if val != '' and val !='\n' and val !='\r\n':
               list_el.append(val)
        for i in range(0,len(list_el)):
             if i!=0:
                 element_ff.append(list_el[i])    
#    print ff,' ',element_ff
        with MPRester(MAPI_KEY) as m:
             data = m.get_entries_in_chemsys(element_ff,inc_structure='final', property_data=["unit_cell_formula","material_id","icsd_id","spacegroup","energy_per_atom","formation_energy_per_atom","pretty_formula","band_gap","total_magnetization","e_above_hull"])
             if (len(element_ff)>1):
                 try:
                     entries = m.get_entries_in_chemsys(element_ff)
                     pd = PhaseDiagram(entries)
                     plotter = PDPlotter(pd, show_unstable=True)
                     image=str(ff)+str("_DFT")+str(".jpg")
                     plotter.write_image(image)
                 except:
                     pass
             structures=[]
             structures_cvn=[]
             icsd_arr=[]
             mp_arr=[]
             sg_arr=[]
             enp_arr=[]
             fenp_arr=[]
             pf_arr=[]
             ucf_arr=[]
             bg_arr=[]
             tm_arr=[]
             ehull_arr=[]
             for d in data:
                 x=d.data["material_id"]
                 sg=d.data["spacegroup"]
                 enp=d.data["energy_per_atom"]
                 fenp=d.data["formation_energy_per_atom"]
                 pf=d.data["pretty_formula"]
                 ucf=d.data["unit_cell_formula"]
                 bg=d.data["band_gap"]
                 tm=d.data["total_magnetization"]
                 ehull=d.data["e_above_hull"]
                 icsd=d.data["icsd_id"]
                 structure = m.get_structure_by_material_id(x)
                 structures.append(structure)
                 icsd_arr.append(icsd)
                 mp_arr.append(x)
                 sg_arr.append(sg)
                 enp_arr.append(enp)
                 fenp_arr.append(fenp)
                 pf_arr.append(pf)
                 bg_arr.append(bg)
                 tm_arr.append(tm)
                 ucf_arr.append(ucf)
                 ehull_arr.append(ehull)
              

                 comment=str("bulk@")+str(x)
                 folder2=str(os.getcwd())+str("/")+str(comment)+str("_fold")
                 if not os.path.exists(folder2):
                     os.makedirs(str(folder2))
                 print ("folder2=",folder2)
                 cwd2=str(os.getcwd())
                 os.chdir(str(folder2))

                 p=Poscar(structure)
                 p.comment=comment
                 p.write_file("POSCAR")
                 poscar_file=str(os.getcwd())+str("/POSCAR")

                 pair_coeff=str(cwd1)+str("/")+str(file)
                 #pair_coeff=str('/data/knc6/JARVIS-FF-NEW/ALLOY')+str("/")+str(file)
                 parameters = {'pair_style':'eam/fs','pair_coeff':pair_coeff,'atom_style': 'charge' ,'control_file':'/users/knc6/inelast.mod'}
                 main_file=open("setup.py","w")
                 line=str("from NEW_LAMMPS10 import main_func")+'\n'
                 main_file.write(line)
                 line=str("from pymatgen.io.vasp.inputs import  Poscar")+'\n'
                 main_file.write(line)
                 #line=str("try:")+'\n'
                 #main_file.write(line)
                 line=str("p=Poscar.from_file(")+str('"')+str(poscar_file)+str('"')+str(")")+'\n'
                 main_file.write(line)
                 line=str("main_func(mat=p")+str(",")+str("parameters=")+str(parameters)+str(")")+'\n'
                 main_file.write(line)
                 main_file.close()
                 cmd=str("nohup python setup.py &")
                 os.system(cmd)
                 #line=str("except:")+'\n'
                 #main_file.write(line)
                 #line=str("     pass")+'\n'
                 #main_file.write(line)
                 #try:
                 #    p=Poscar.from_file(poscar_file)
                 #    main_func(mat=p,parameters=parameters)
                 #except:
                 #     pass
                 os.chdir(cwd2)#=str(os.getcwd())
        os.chdir(cwd1)#=str(os.getcwd())
    except:
       pass
