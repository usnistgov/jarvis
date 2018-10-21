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
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from  jarvis.lammps.jlammps import vac_antisite_def_struct_gen,surfer
import numpy as np,time,json
import sys,os,subprocess,socket
from pymatgen.io.ase import AseAtomsAdaptor
#from ase.calculators.lammpsrun import LAMMPS, Prism
import sys,zipfile
import fortranformat as fform
from pymatgen.core.structure import Structure
import os
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
import sys,shutil
#import urllib2
#import urllib
#from BeautifulSoup import BeautifulSoup
import json
MAPI_KEY = os.environ.get("MAPI_KEY", "")
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from math import ceil
from collections import OrderedDict

import matplotlib
matplotlib.use('Agg')

from pymatgen.core.composition import Composition
#from pymatgen.phasediagram.entries import PDEntry
#from pymatgen.phasediagram.pdmaker import PhaseDiagram
#from pymatgen.phasediagram.plotter import PDPlotter

import glob
try:
  from ase.calculators.lammpsrun import LAMMPS, Prism
except:
   print ('Please install ase')
   pass

###################MAIN PART############################################
for file in glob.glob("CdTeZnSeHgS0.sw"):
#for file in glob.glob("ff.MoS.10.26.2016.comb"):
    try:
        folder1=str(os.getcwd())+str("/")+str(file)+str("_nist")
        if not os.path.exists(folder1):
           os.makedirs(str(folder1))
        cwd1=str(os.getcwd())
        print ("folder1=",folder1)
        ff=str(file)
        element_ff=[]
        #f=open(ff,'r')
        os.chdir(str(folder1))
        #list_el=[]
        #lines=f.readlines()
        #content=(lines[3]).split(" ")
        #content=(lines[3]).split("' '|\n|\r\n")
        #for val in content:
        
       #     if val != '' and val !='\n' and val !='\r\n':
       #        list_el.append(val)
       # for i in range(0,len(list_el)):
       #      if i!=0:
       #          element_ff.append(list_el[i])    
       #    print ff,' ',element_ff
        element_ff=['Cd','Te','Zn','Se','Hg','S']
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
                 parameters = {'job_bin':'/cluster/bin/lmp_ctcms-14439-knc6-2<in.elastic>outt','pair_style':'sw','pair_coeff':pair_coeff,'atom_style': 'charge' ,'control_file':'/users/knc6/inelast.mod'}

                 f=open('setup.py','w')
                 line=str("from NEW_LAMMPS10a import main_func")+'\n'
                 f.write(line)
                 line=str("from pymatgen.io.vasp.inputs import  Poscar")+'\n'
                 f.write(line)
                 line=str("p=Poscar.from_file(")+str('"')+str(poscar_file)+str('"')+str(")")+'\n'
                 f.write(line)
                 line=str("main_func(mat=p")+str(",")+str("parameters=")+str(parameters)+str(")")+'\n'
                 f.write(line)
                 f.close()

                 nprocs=1
                 nnodes=1
                 
                 f=open("submit_job","w")
                 line=str("#!/bin/bash")+'\n'
                 f.write(line)
                 line=str("#PBS -N ")+str(''.join(element_ff))+str(x)+'\n'
                 f.write(line)
                 line=str("#PBS -o test.log")+'\n'
                 f.write(line)
                 line=str("#PBS -m abe")+'\n'
                 #f.write(line)
                 line=str("#PBS -j oe")+'\n'
                 f.write(line)
                 line=str("#PBS -q highmem")+'\n'
                 f.write(line)
                 line=str("#PBS -r n")+'\n'
                 f.write(line)
                 line=str("#PBS -l nodes=")+str(nnodes)+str(":")+str("ppn=")+str(int(float(nprocs)/float(nnodes)))+'\n'
                 f.write(line)
                 dir=str(os.getcwd())
                 line=str("cd ")+dir+'\n'
                 f.write(line)
                 line=str("/users/knc6/anaconda2/bin/python setup.py>sout")+'\n'
                 f.write(line)
                 f.close()
 

                 with open('job.out', 'w') as f:
                    p = subprocess.Popen(['qsub','submit_job'], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    stdout, stderr = p.communicate()
                    job_id = stdout.rstrip('\n').split()[-1]
                    print ("stdout,stderr",stdout, stderr)
                    #job_id = str(stdout.split('Your job')[1].split(' ')[1])
                    f.write(job_id)

                 os.chdir(cwd2)#=str(os.getcwd())
        os.chdir(cwd1)#=str(os.getcwd())
    except:
       pass
