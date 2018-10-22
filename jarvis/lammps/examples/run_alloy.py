import glob
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import  Poscar
from pymatgen.matproj.rest import MPRester
from pymatgen.core.structure import Structure
import os
import json
MAPI_KEY = os.environ.get("MAPI_KEY", "")
from pymatgen.matproj.rest import MPRester
import matplotlib,subprocess
matplotlib.use('Agg')

from pymatgen.core.composition import Composition

def your_job(name=''):
        nprocs=1
        nnodes=1

        f=open('submit_job','w')
        line=str("#!/bin/bash")+'\n'
        f.write(line)
        line=str("#SBATCH -J ")+str(name)+'\n'
        f.write(line)
        #line=str("#SBATCH -t 330:10:00")+'\n'
        #f.write(line)
        line=str("#SBATCH -o test.log")+'\n'
        f.write(line)
        line=str("#SBATCH -N ")+str(nnodes)+'\n'
        #f.write(line)
        line=str("#SBATCH --ntasks-per-node=")+str(nprocs)+'\n'
        f.write(line)
        #line=str("#SBATCH -p mml")+'\n'
        #f.write(line)
        #line=str("#SBATCH --mem=")+str(mem)+'\n'
        #f.write(line)
        dir=str(os.getcwd())
        line=str("cd ")+str(dir)+'\n'
        f.write(line)
        line=str(" source ~/anaconda2/envs/my_jarvis/bin/activate my_jarvis")+'\n'
        f.write(line)
        line=str('python setup.py>out')+'\n'
        f.write(line)
        f.close()
        with open('jff.out', 'w') as f:
                    p = subprocess.Popen(['sbatch','submit_job'], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    stdout, stderr = p.communicate()
                    job_id = str(stdout.split('.')[0])
                    print ("stdout,stderr",stdout, stderr)
                    f.write(job_id)
   


for file in glob.glob("*.alloy"):
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
        for val in content:
        
            if val != '' and val !='\n' and val !='\r\n':
               list_el.append(val)
        for i in range(0,len(list_el)):
             if i!=0:
                 element_ff.append(list_el[i])    
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
             for d in data:
                 x=d.data["material_id"]
                 structure = m.get_structure_by_material_id(x)
              

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
                 parameters = {'phonon_control_file':'/users/knc6/in.phonon','surf_control_file':'/users/knc6/inelast_nobox.mod','def_control_file':'/users/knc6/inelast_nobox.mod','json_dat':'/rk2/knc6/JARVIS-DFT/MDCS/all_mp.json','c_size':3,'vac_size':25,'surf_size':25,'phon_size':20,'cluster':'sbatch','exec':'/users/knc6/Software/LAMMPS/lammps-master/src/lmp_serial <in.main >out','pair_style':'eam/alloy','pair_coeff':pair_coeff,'atom_style': 'charge' ,'control_file':'/users/knc6/inelast.mod'}
                 main_file=open("setup.py","w")
                 line=str("from jarvis.lammps.jlammps import main_func")+'\n'
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
                 try:
                   your_job(name=x)
                 except:
                  pass
                 os.chdir(cwd2)
                 
        os.chdir(cwd1)
    except:
       pass
#your_job('ll')
