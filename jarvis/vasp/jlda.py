"""
Module to run LDA based High-throughput calculations
"""
from __future__ import division, unicode_literals, print_function
import os
from monty.json import MontyEncoder, MontyDecoder
from custodian.vasp.jobs import VaspJob
from pymatgen.io.vasp import VaspInput, Vasprun
#from mpinterfaces.data_processor import MPINTVasprun
#from mpinterfaces.instrument import MPINTVaspInputSet, MPINTVaspJob
from custodian.vasp.jobs import VaspJob
from pymatgen.io.vasp.outputs import Oszicar
from pymatgen.io.vasp.outputs import Vasprun
from subprocess import Popen, PIPE
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import sys,shutil,glob,codecs
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import  Slab, SlabGenerator, generate_all_slabs,get_symmetrically_distinct_miller_indices
from custodian.vasp.handlers import VaspErrorHandler, UnconvergedErrorHandler, \
    MeshSymmetryErrorHandler, NonConvergingErrorHandler, PotimErrorHandler
import json,yaml
from numpy import linalg as LA
import time
from collections import OrderedDict
#from  gen_def_surf import vac_antisite_def_struct_gen,surfer
from pymatgen.ext.matproj import MPRester
import subprocess
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar, Poscar, VaspInput
from pymatgen.io.vasp.inputs import Potcar, Kpoints
#from pymatgen.io.vaspio_set import MPVaspInputSet, MPNonSCFVaspInputSet
from numpy import matrix
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from ase.lattice.cubic import DiamondFactory, SimpleCubicFactory
from ase.lattice import bulk
from ase.lattice.compounds import Zincblende,AuCu3,Rocksalt,AuCu,CsCl,TRI_Fe2O3,HEX_Fe2O3,CsClFactory
from ase.lattice.spacegroup import crystal
import operator



def  make_prototypes(strt=None):
     strt.sort()
     symbs=[str(specie.symbol) for specie in strt.types_of_specie]
     protos=[]
     if len(symbs)==1:
     # FCC, BCC,SC,HCP,DIAM
        hcp = bulk(symbs, 'hcp',a=4.0)
        hcp1= AseAtomsAdaptor().get_structure(hcp)
        sg_mat = SpacegroupAnalyzer(hcp1)
        hcp2 = Poscar(sg_mat.get_conventional_standard_structure())
        hcp2.comment=("bulk@HCP")
        print ("HCP",hcp2)
        protos.append(hcp2)
        sc = bulk(symbs, 'sc', a=2.8)
        sc1= AseAtomsAdaptor().get_structure(sc)
        sg_mat = SpacegroupAnalyzer(sc1)
        sc2 = Poscar(sg_mat.get_conventional_standard_structure())
        sc2.comment=("bulk@SC")
        print ("SC",sc2)
        protos.append(sc2)
        fcc = bulk(symbs, 'fcc', a=4.0)
        fcc1= AseAtomsAdaptor().get_structure(fcc)
        sg_mat = SpacegroupAnalyzer(fcc1)
        fcc2 = Poscar(sg_mat.get_conventional_standard_structure())
        fcc2.comment=("bulk@FCC")
        print ("FCC",fcc2)
        protos.append(fcc2)
        bcc = bulk(symbs, 'bcc', a=3.3)
        bcc1= AseAtomsAdaptor().get_structure(bcc)
        sg_mat = SpacegroupAnalyzer(bcc1)
        bcc2 = Poscar(sg_mat.get_conventional_standard_structure())
        bcc2.comment=("bulk@BCC")
        print ("BCC",bcc2)
        protos.append(bcc2)



        diam = bulk(symbs[0], 'diamond',a=6.0)
        diam1= AseAtomsAdaptor().get_structure(diam)
        sg_mat = SpacegroupAnalyzer(diam1)
        diam2 = Poscar(sg_mat.get_conventional_standard_structure())
        diam2.comment=("bulk@DIAM")
        print ("DIAM",diam2)
        protos.append(diam2)

     elif len(symbs)==2:
     # RS,ZB,CsCl,TiO2,CaF2
        rc=Rocksalt(directions=[[1,0,0], [0,1,0], [0,0,1]],
            size=(1,1,1),
            symbol=symbs,
            latticeconstant=4.50)
        rc1= AseAtomsAdaptor().get_structure(rc)
        sg_mat = SpacegroupAnalyzer(rc1)
        rc2 = Poscar(sg_mat.get_conventional_standard_structure())
        rc2.comment=("bulk@RS")
        print (rc2)
        protos.append(rc2)

        zb=Zincblende(directions=[[1,0,0], [0,1,0], [0,0,1]],
            size=(1,1,1),
            symbol=symbs,
            latticeconstant=4.50)
        zb1= AseAtomsAdaptor().get_structure(zb)
        sg_mat = SpacegroupAnalyzer(zb1)
        zb2 = Poscar(sg_mat.get_conventional_standard_structure())
        zb2.comment=("bulk@ZB")
        print (zb2)
        protos.append(zb2)

        cscl=CsCl(directions=[[1,0,0], [0,1,0], [0,0,1]],
            size=(1,1,1),
            symbol=symbs,
            latticeconstant=4.50)
        cscl1= AseAtomsAdaptor().get_structure(cscl)
        sg_mat = SpacegroupAnalyzer(cscl1)
        cscl2 = Poscar(sg_mat.get_conventional_standard_structure())
        cscl2.comment=("bulk@CsCl")
        print (cscl2)
        protos.append(cscl2)

        rutile =crystal(symbs, basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                spacegroup=136, cellpar=[4.6, 4.6, 2.95, 90, 90, 90])
        rut1= AseAtomsAdaptor().get_structure(rutile)
        sg_mat = SpacegroupAnalyzer(rut1)
        rut2 = Poscar(sg_mat.get_conventional_standard_structure())
        rut2.comment=("bulk@Rutile")
        print (rut2)
        protos.append(rut2)

        joined=str(''.join(symbs))+str(2)
        print (joined)
        caf2 = bulk((joined), 'fluorite', a=4.0)
        caf21= AseAtomsAdaptor().get_structure(caf2)
        sg_mat = SpacegroupAnalyzer(caf21)
        caf22 = Poscar(sg_mat.get_conventional_standard_structure())
        rut2.comment=("bulk@Fluorite")
        print (caf22)
        protos.append(caf22)
     else:
        print ("Add case")
     return protos

def smart_protos(protos=[]):
    for p in protos:
           enp,contc=smart_converge(mat=p)
           #pass
def sub_job(file=None):
    process = Popen(['qsub','-cwd', '-pe', 'nodal', '16', file], stderr=PIPE)

      
def check_polar(file):
    up=0
    dn=0
    coords = np.array(file.frac_coords)
    z_max = max(coords[:, 2])
    z_min=  min(coords[:, 2])
    for site in file:
        if (site.frac_coords[2]==z_max) :
           up=up+site.specie.number
        if (site.frac_coords[2]==z_min) :
           dn=dn+site.specie.number
    polar=False
    if(up !=dn):
       print ('polar')
       polar=True
    if(up == dn):
       print ('Non-polar')
       polar=False
    return polar


def getvasp(strt=None,encut=600,kpoints=None):
#def getvasp(strt):
    #Get All the INCAR,KPOINTS,POTCAR
    v = MPVaspInputSet(force_gamma='True')
    #strt=Structure.from_file('POSCAR')
    v.write_input(strt,'./')
    incar_dict = { 'SYSTEM': 'slab',
                   'ENCUT': encut,
                   'AlGO': 'accurate',
                   'ISMEAR': 0,
                   'ISYM':0,
                   'ADDGRID': '.TRUE.',
                   'EDIFF': 1e-09,
                   'NPAR': 4,
                   'SIGMA': 0.1,
                   'LORBIT': 11,
                   'PREC': 'Accurate',
                   'LWAVE':  '.FALSE.',
                   'LCHARG': '.FALSE.',
               }
    incar = Incar.from_dict(incar_dict)
    incar.write_file("INCAR")
    #num_kpoints = den*strt.lattice.reciprocal_lattice.volume
    #kpoints = Kpoints.automatic_density(strt, num_kpoints * strt.num_sites,force_gamma=True)
    kpoints.write_file('KPOINTS')
def pos2phonts(strt=None,mesh=[35,35,35]):
    #Getting conventional cell#
    #strt = Structure.from_file("POSCAR")
    #strt=strt.structure
    sg_mat = SpacegroupAnalyzer(strt)
    strt = sg_mat.get_conventional_standard_structure()
    f=open('phonons_input.dat','w')
    #line=str('kpoints 35 35 35 1')+'\n'
    line=str('kpoints')+str(' ')+str(mesh[0])+str(' ')+str(mesh[1])+str(' ')+str(mesh[2])+str(' 1')+'\n'
    f.write(line)
    line=str('phonons_only T')+'\n'
    f.write(line)
    line=str('pdos 0 20 300 10')+'\n'
    f.write(line)
    #Just to look separated
    line='\n'
    f.write(line)
    line='\n'
    f.write(line)
    line='\n'
    f.write(line)
    line='\n'
    f.write(line)
    line=str('vectors')+'\n'
    f.write(line)
    line=str(round(float(strt.lattice.matrix[0][0]),4))+'  '+str(round(float(strt.lattice.matrix[0][1]),4))+'  '+str(round(float(strt.lattice.matrix[0][2]),4))+'\n'
    f.write(line)
    line=str(round(float(strt.lattice.matrix[1][0]),4))+'  '+str(round(float(strt.lattice.matrix[1][1]),4))+'  '+str(round(float(strt.lattice.matrix[1][2]),4))+'\n'
    f.write(line)
    line=str(round(float(strt.lattice.matrix[2][0]),4))+'  '+str(round(float(strt.lattice.matrix[2][1]),4))+'  '+str(round(float(strt.lattice.matrix[2][2]),4))+'\n'
    f.write(line)
#Just to look separated
    line='\n'
    f.write(line)
    line='\n'
    f.write(line)
    line='\n'
    f.write(line)
    line='\n'
    f.write(line)
    line=str('species')+'  '+str(len(strt.composition.elements))+'\n'
    f.write(line)
    for i in  range (len(strt.composition.elements)):
        mass=str(strt.composition.elements[i].atomic_mass)
        mass = mass.replace('amu', ' 0.0')
        line=str(strt.composition.elements[i])+'  '+mass+'\n'
        f.write(line)
    line=str('Lattice 1.0')+'\n'
    f.write(line)

    line='\n'
    f.write(line)
    line='\n'
    f.write(line)
    line=str('AbInitio T F')+'\n'
    f.write(line)
    line=str('FP_interface VASP')+'\n'
    f.write(line)
    line=str('D3_cutoff 4.0')+'\n'
    f.write(line)
    line=str('delta 0.05')+'\n'
    f.write(line)
    line=str('numerical_2der T')+'\n'
    f.write(line)
    line=str('natoms')+' '+str(int(strt.composition._natoms))+'\n'
    f.write(line)
    line=str('fractional')+'\n'
    f.write(line)


    for site in strt:
#    print site.frac_coords,site.specie,site.specie.number
        line= str(site.specie)+'  '+str(site.frac_coords[0])+'  '+str(site.frac_coords[1])+'  '+str(site.frac_coords[2])+'\n'
        f.write(line)

    line=str('end')+'\n'
    f.write(line)


def Raman(strt=None,encut=500,length=50,efield='0.0 0.0 0.0'):
    kpoints=Auto_Kpoints(mat=strt,length=length)
    RAMANDIR=str("RAMANDIR-")+str(strt.comment)
    comm=str(strt.comment)
    strt=strt.structure
    if not os.path.exists(RAMANDIR):
       os.makedirs(RAMANDIR)
    os.chdir(RAMANDIR)
    strt.sort()


   # v = MPVaspInputSet(force_gamma='True')
   # v.write_input(s,DOSDIR)
    incar_dict = { 'SYSTEM': 'PHONON',
                   'ENCUT': encut,
                   'AlGO': 'Normal',
                   'IBRION': 8,
                   'POTIM': 0.01,
                   'LPLANE': '.FALSE.',
                   'LWAVE': '.FALSE.',
                   'LCHARG': '.FALSE.',
                   'LELF': '.FALSE.',
                   'LVTOT': '.FALSE.',
                   'ISTART': 0,


                   'NWRITE': 3,
                   'ISMEAR': 0,
                   'EDIFF': 1e-08,
                   'NPAR': 16,
                   'SIGMA': 0.05,
                   'LREAL': '.FALSE.',
                   'ADDGRID': '.TRUE.',
                   'LEPSILON': '.TRUE.',
                   'LORBIT': 11,
                   'PREC': 'Accurate',
               }
    incar = Incar.from_dict(incar_dict)
    #incar.write_file("INCAR")
    #kpoints=Kpoints.from_file('/scratch/lfs/kamal/VASP/Raman/ZnO-Ang/ZnO222/KPOINTS')
    #num_kpoints = den*strt.lattice.reciprocal_lattice.volume
    #kpoints = Kpoints.automatic_density(strt, num_kpoints * strt.num_sites,force_gamma=True)
    #kpoints.write_file('KPOINTS')
    c_size=12
    r1=int(round(float(c_size)/float( max(abs(strt.lattice.matrix[0])))))
    r2=int(round(float(c_size)/float( max(abs(strt.lattice.matrix[1])))))
    r3=int(round(float(c_size)/float( max(abs(strt.lattice.matrix[2])))))
    if r1<1:
       r1=1
    if r2<1:
       r2=1
    if r3<1:
       r3=1
    strt.make_supercell([r1,r2,r3])
#    strt.make_supercell([dim1,dim2,dim3])

#    if comm.startswith('Surf') and float( max(strt.lattice.matrix[0])) <11.0 and float( max(strt.lattice.matrix[1])) <11.0 :
#       c_size=12
 #      dim1=int((float(c_size)/float( max(strt.lattice.matrix[0]))))
 #      dim2=int(float(c_size)/float( max(strt.lattice.matrix[1])))
 #      strt.make_supercell([dim1,dim2,1])
 #   elif float( max(strt.lattice.matrix[0])) <11.0 and float( max(strt.lattice.matrix[1])) <11.0 and float( max(strt.lattice.matrix[2])) <11.0:
 #      c_size=12
 #      dim1=int((float(c_size)/float( max(strt.lattice.matrix[0]))))
 #      dim2=int(float(c_size)/float( max(strt.lattice.matrix[1])))
 #      dim3=int(float(c_size)/float( max(strt.lattice.matrix[2])))
 #      strt.make_supercell([dim1,dim2,dim3])

    en,contc=run_job(mat=Poscar(strt),incar=incar,kpoints=kpoints,jobname=str('PHONON'))
    path=str(contc.split('/CONTCAR')[0])
    print ("copying files")
    cmd1=str('rm -f KPOINTS POTCAR OUTCAR OUTCAR.phon POSCAR POSCAR.phon')
    os.system(cmd1)
    cmd1=str('ln -s ')+str(path)+'/KPOINTS ./KPOINTS'
    os.system(cmd1)
    cmd1=str('ln -s ')+str(path)+'/POTCAR ./POTCAR'
    os.system(cmd1)
    cmd1=str('ln -s ')+str(path)+'/POSCAR ./POSCAR.phon'
    os.system(cmd1)
    cmd1=str('ln -s ')+str(path)+'/OUTCAR ./OUTCAR.phon'
    os.system(cmd1)

        #user_incar_settings={"EDIFF":1E-4,"NBANDS":30,"ISIF":2,"NSW":0}
    #mpvis = MPNonSCFVaspInputSet(user_incar_settings=user_incar_settings)
    os.environ['VASP_RAMAN_PARAMS'] = '01_18_2_0.01'
    os.environ['VASP_RAMAN_RUN'] = 'mpirun -np 16 /users/knc6/VASP/vasp54/src/vasp.5.4.1/bin/vasp_std >vasp.out'
    #mpvis.write_input(s,BANDSDIR)
    if efield!='0.0 0.0 0.0':
       incar_dict = { 'SYSTEM': 'PHONON',
                   'ENCUT': encut,
                   'AlGO': 'Normal',
                   'ISYM': 0,
                   'LPLANE': '.FALSE.',
                   'LWAVE': '.FALSE.',
                   'LCHARG': '.FALSE.',
                   'LELF': '.FALSE.',
                   'LVTOT': '.FALSE.',
                   'ISTART': 0,
                   'NWRITE': 3,
                   'ISMEAR': 0,
                   'LEPSILON':'.TRUE.',
                   'LCALCEPS':'.TRUE.',
                   'EFIELD_PEAD': efield,

                   'EDIFF': 1e-08,
                   'KPAR': 4,
                   'SIGMA': 0.05,
                   'LREAL': '.FALSE.',
                   'ADDGRID': '.TRUE.',
                   'LORBIT': 11,
                   'PREC': 'Accurate',
               }
    else:
       incar_dict = { 'SYSTEM': 'PHONON',
                   'ENCUT': encut,
                   'AlGO': 'Normal',
                   'ISYM': 0,
                   'LPLANE': '.FALSE.',
                   'LWAVE': '.FALSE.',
                   'LCHARG': '.FALSE.',
                   'LELF': '.FALSE.',
                   'LVTOT': '.FALSE.',
                   'ISTART': 0,
                   'NWRITE': 3,
                   'ISMEAR': 0,
                   'LEPSILON':'.TRUE.',

                   'EDIFF': 1e-08,
                   'KPAR': 4,
                   'SIGMA': 0.05,
                   'LREAL': '.FALSE.',
                   'ADDGRID': '.TRUE.',
                   'LORBIT': 11,
                   'PREC': 'Accurate',
               }
 
    incar = Incar.from_dict(incar_dict)
    incar.write_file("INCAR")
    kpoints.write_file('KPOINTS')


    f=open('raman.job','w')
    try:
       line=str("#PBS -N Phonon-")+str((Poscar(strt)).comment)+'\n'
    except:
       line=str("#PBS -N Phonon")+'\n'
    f.write(line)
    line=str("#PBS -m abe")+'\n'
    f.write(line)
    line=str("#PBS -W group_list=ssdarpa")+'\n'
    f.write(line)
    line=str("#PBS -j oe")+'\n'
    f.write(line)
    line=str("#PBS -r n")+'\n'
    f.write(line)
    line=str("#PBS -l nodes=4:ppn=8")+'\n'
    f.write(line)
    line=str("#PBS -l pmem=1500mb")+'\n'
    f.write(line)
    line=str("#PBS -o test.log")+'\n'
    f.write(line)
    line=str("#PBS -l walltime=48:10:00")+'\n'
    f.write(line)
    dir=str(os.getcwd())
    line=str("cd ")+dir+'\n'
    f.write(line)
    line=str("export VASP_RAMAN_RUN='mpirun -np 16 /users/knc6/VASP/vasp54/src/vasp.5.4.1/bin/vasp_std >vasp.out'")+'\n'
    f.write(line)
    nmodes=3*int(strt.num_sites)
    input=str("'")+str("01_")+str(nmodes)+str("_2_0.01'")
    line=str("export VASP_RAMAN_PARAMS=")+input+'\n'
    #line=str("export VASP_RAMAN_PARAMS='01_12_2_0.01'")+'\n'
    f.write(line)
    line=str("python /users/knc6/bin/vasp_raman.py > vasp_raman.out")+'\n'
    f.write(line)
    
    f.close()

    #shutil.copy2('/data/knc6/Arunima/Arunima/1T_prime/vdw_kernel.bindat','./')
    sub_job(file='raman.job')
#    strt=get_cvn(strt)
#    pos2phonts(strt=strt)
    os.chdir('../')
def _submit_to_queue(script_file):
        # submit the job, return process and pid.
        process = Popen(("/bin/bash", script_file), stderr=PIPE)
        return SubmitResults(qid=process.pid, out='no out in shell submission', err='no err in shell submission', process=process)
    

def phonts(strt=None,encut=600,kpoints=None):
    comm=str(strt.comment)
    PHONTSDIR=str("PHONTSDIR-")+str(comm)
    main_json=str("MAIN-RELAX-")+str(comm)
   # for file in glob.glob("*.json"):
   #     if file.startswith(main_json):
   #       den =int(file.split("-")[-1])
   #       main_file=jobs_from_file(file)
   #       for job in main_file:
   #           encut=float(job.vis.incar["ENCUT"])
    strt=strt.structure
    os.makedirs(PHONTSDIR)
    #strt=Structure.from_file('POSCAR')
    strt.sort()

    os.chdir(PHONTSDIR)
#    strt=get_cvn(strt)
    pos2phonts(strt=strt)
    shutil.copy2('phonons_input.dat','phonons_input.dat_1')
    exe = str("PhonTS phonons_input.dat")
    os.system(exe)
    for b in glob.glob('*.str'):
        #cnt=cnt+1
        folder=str(b+'.dir')
        print (folder)
        #if not os.path.exists(folder):
        os.makedirs(folder)
        shutil.copy2(b,folder)
        os.chdir(folder)
        print ("folder===",folder)
        shutil.copy2('/users/knc6/phonon.job','phonon.job')
        super=Poscar.from_file(b)
        getvasp(super.structure,kpoints=kpoints,encut=encut)
        shutil.copy2(b,'POSCAR')
        print ("bbbbbbbbbbbbbbbbbbbbbbb=",b)
        exe = str("qsub -cwd -pe nodal 16 phonon.job")
        #exe = str("mpiexec vasp")
        os.system(exe)
        filename = b
        (prefix, sep,suff) = filename.rpartition('.str')
        new_filename = prefix + '.out'
        status=False
        while status !=True:
            try:
                #file=os.path.join(str(filename)+str('.dir'),'/vasprun.xml')
                file=str(os.getcwd())+str('/vasprun.xml')
                print (file)
                run=Vasprun(file)
                print ("file=",run)
                status=run.converged
                print ("status=",status)
                time.sleep(5)
            except:
                time.sleep(5)
        inp=str('../')+str(new_filename)
        shutil.copy2('OUTCAR',inp)
        os.chdir('../')
    replceTF()
    exe = str("PhonTS phonons_input.dat")
    os.system(exe)
    os.chdir('../')

def replceTF():
    f = codecs.open('phonons_input.dat',encoding='utf-8')
    contents = f.read()
    newcontents = contents.replace('AbInitio T F','AbInitio F T')
    f1=open('phonons_input.dat','w')
    f1.write(newcontents)
    f1.close()
def get_lowest_en_from_mp(formula, MAPI_KEY="", all_structs=False):
    if not MAPI_KEY:
        MAPI_KEY = os.environ.get("MAPI_KEY", "")
        if not MAPI_KEY:
            print('API key not provided')
            print('get API KEY from materialsproject and set it to the MAPI_KEY environment variable. aborting ... ')
            sys.exit()
    with MPRester(MAPI_KEY) as m:
        data = m.get_data(formula)
        structures = []
        x = {}
        print("\nnumber of structures matching the chemical formula {0} = {1}".format(formula, len(data)))
        print("The one with the the lowest energy above the hull is returned, unless all_structs is set to True")
        for d in data:
            mpid = str(d['material_id'])
            x[mpid] = d['e_above_hull']
            if all_structs:
                structure = m.get_structure_by_material_id(mpid)
                structures.append(structure)
        else:
            mineah_key = sorted(x.items(), key=operator.itemgetter(1))[0][0]
            print("The id of the material corresponding to the lowest energy above the hull = {0}".format(mineah_key))
            if mineah_key:
                with MPRester(MAPI_KEY) as m:
                    data = m.get_data(mineah_key)
                    x = {}
                    for d in data:
                         x['energy_per_atom'] = str(d['energy_per_atom'])
                         enp= x['energy_per_atom']
                #return m.get_structure_by_material_id(mineah_key)
                return enp
            else:
                return None

def sum_chem_pot(strt=None):
    sum=0
    symb=strt.symbol_set
    for el in symb:
        enp=get_lowest_en_from_mp(el)
        sum=float(sum)+float(enp)
    return sum
#strt=Structure.from_file('/home/kamal/POSCAR')
#sum=sum_chem_pot(strt)
#print sum

# all the info/warnings/outputs redirected to the log file: convg.log

poscar_list = []

def run_cal(turn_knobs, qadapter, job_cmd, job_dir, checkpoint_file,
            incar=None, poscar=None, potcar=None, kpoints=None,
            Grid_type='G',functional='LDA',is_matrix=True):
    handlers = []
    outfile=str(os.getcwd())+str('/')+str('vasp.out')
    handlers = [VaspErrorHandler(output_filename=outfile)] #, MeshSymmetryErrorHandler(),
    #handlers = [VaspErrorHandler(), MeshSymmetryErrorHandler(),
    #            UnconvergedErrorHandler(), NonConvergingErrorHandler(),
    #            PotimErrorHandler()]
    cal = Calibrate(incar, poscar, potcar, kpoints,
                    is_matrix=is_matrix, 
                    turn_knobs=turn_knobs,handlers=handlers, qadapter=qadapter,
                    job_cmd = job_cmd, job_dir=job_dir,
                    Grid_type=Grid_type,functional=functional,
                    checkpoint_file=checkpoint_file, cal_logger=logger)
    cal.setup()
    #c = Custodian(handlers, cal, max_errors=5)
    #print ("custodian",c,type(c))
    #c.run()
    cal.run()
def run_job(mat=None,incar=None,kpoints=None,jobname=''):    
     
    poscar_list=[(mat)]
    nprocs = 16
    nnodes = 8
    mem='1500'
    job_name=str(mat.comment)
    walltime = '148:12:00'
    #job_bin = 'mpirun -np 16 /users/knc6/VASP/vasp54/src/vasp.5.4.1/bin/vasp_std >vasp.out'
    job_bin =  'mpirun -np 16 /users/knc6/bin/vasp_5.3.3_parallel >vasp.out'
    #qadapter, job_cmd = get_run_cmmnd(nnodes=nnodes, nprocs=nprocs,
    #                                  walltime=walltime,
    #                                  job_bin=job_bin, mem=mem,job_name=job_name)
    #turn_knobs = OrderedDict([('POSCAR', poscar_list)])
    job_dir=str(jobname)
    checkpoint_files = []
    chkpt_file = str(jobname)+str('.json')
    checkpoint_files.append(chkpt_file)
    if mat.comment.startswith('Surf'):
       job_bin =  'mpirun -np 16 /data/knc6/Arunima/Arunima-Kamal/vasp_no_z_vdW  >vasp.out'
       #job_bin =  'mpirun -np 16 /users/knc6/VASP/vasp54/src/vasp.5.4.1noz/bin/vasp_std >vasp.out'
       [a,b,c]=kpoints.kpts[0]
       kpoints.kpts=[[a,b,1]]
       try:
           pol=check_polar(mat.structure)
           if pol==True:
                ase_atoms = AseAtomsAdaptor().get_atoms(mat.sttucture)
                COM=ase_atoms.get_center_of_mass(scaled=True)
                print ("COM=",COM)
                incar.update({"LDIPOL": '.TRUE.',"IDIPOL":3,"ISYM": 0,"DIPOL":COM})
                print ("Polar surface encountered in run_job",mat.comment)
       except:
            pass
    with open('/users/knc6/bin/Special_POTCAR.yaml', 'r') as f:
         doc = yaml.load(f)
         pots=doc['POTCAR']
    new_symb=[]
    for el in mat.site_symbols:
        new_symb.append(pots[el])
    potcar = Potcar(symbols=new_symb,functional="LDA")
    wait=False
    if not os.path.exists(jobname):
      need_run = True
      os.makedirs(jobname)
      os.chdir(jobname)
      #shutil.copy2('/data/knc6/Arunima/Arunima/1T_prime/vdw_kernel.bindat','./')
      f=open('submit_job','w')
      line=str("#!/bin/bash")+'\n'
      f.write(line)
      line=str("#PBS -N ")+jobname+'\n'
      f.write(line)
      line=str("#PBS -l walltime=148:10:00")+'\n'
      f.write(line)
      line=str("#PBS -o test.log")+'\n'
      f.write(line)
      line=str("#PBS -m abe")+'\n'
      f.write(line)
      line=str("#PBS -j oe")+'\n'
      f.write(line)
      line=str("#PBS -r n")+'\n'
      f.write(line)
      line=str("#PBS -l nodes=")+str(nnodes)+str(":")+str("ppn=")+str(int(float(nprocs)/float(nnodes)))+'\n'
      f.write(line)
      line=str("#PBS -l pmem=")+str(mem)+'\n'
      f.write(line)
      dir=str(os.getcwd())
      line=str("cd ")+dir+'\n'
      f.write(line)
      line=str(job_bin)+'\n'
      f.write(line)
      f.close()
      incar.write_file("INCAR")
      potcar.write_file("POTCAR")
      kpoints.write_file("KPOINTS")
      mat.write_file("POSCAR")
      import sys
      vinput = VaspInput.from_directory(".")
      job=VaspJob(['qsub -cwd -pe nodal 16','submit_job'], final=False, backup=False)
    #handlers = [VaspErrorHandler(), UnconvergedErrorHandler()]
      handlers = [VaspErrorHandler(), MeshSymmetryErrorHandler(),
		UnconvergedErrorHandler(), NonConvergingErrorHandler(),
		PotimErrorHandler()]
      ticks=time.time()
      print ("job start time for ", jobname,"=",ticks)
    #c = Custodian(handlers, [job],max_errors=1)
    #c.run()
      with open('job.out', 'w') as f:
                p = subprocess.Popen(['qsub','-cwd', '-pe', 'nodal', '16', 'submit_job'], stdout=subprocess.PIPE,
				     stderr=subprocess.PIPE)
                stdout, stderr = p.communicate()
		#self.job_id = stdout.rstrip('\n').split()[-1]
                print ("stdout,stderr",stdout, stderr)
                job_id = str(stdout.split('Your job')[1].split(' ')[1])
                f.write(job_id)
    #output = subprocess.check_output(['qstat', '-j', job_id])
      need_subrun=True
      while need_subrun !=False:
         try:
            output=subprocess.check_output([str('qstat'), str('-j'), job_id])
            time.sleep(5)
         except:
             need_subrun=False
             pass
       
    else:
      os.chdir(jobname)
      print ("directory is CUTION",os.getcwd())
      wait=False
      try:
         jobout=open('job.out','r')
         for jline in jobout:
             jid=(str(jline)).split('\n')[0]
         output=subprocess.check_output([str('qstat'), str('-j'), jid])
         print ('jobout' , jid)
         import sys
         print ("exiting")
         wait=True
         #sys.exit()
         #raise SystemExit
      except:
           need_run=False
           pass
      print ("need_run is now",need_run)
    #correc=open("corrections",'w')
    
    print ("GOES HERE111111111111111111")
    try:
    #if os.path.isfile("OSZICAR"):
    #    oszicar = Oszicar("OSZICAR")
    #    ionic_steps=len(oszicar.ionic_steps)
       if os.path.isfile("OUTCAR"):
    
          wait=Vasprun("vasprun.xml").converged
       #f=open("OUTCAR",'r')
       #for line in f:
       #    if "NSW" in line:
       #        nsw= int(str(line.split("=")[1]).split("number of steps")[0])
       #print ("Going here1",nsw,ionic_steps)
       #for line in f:
       # 	   if "General timing and accounting informations for this job" in line :
       #                 print (line)
       # 	        if  int(nsw)>=int(ionic_steps):
       #                    print ("WAIT=")
       #                    wait=True
       #f.close()
    except:
          pass
    msg=[] #Error messagesMPInterfaces: A Materials Project based Python Tool for High-Throughput Computational Screening of Interfacial Systems
    kp_old=Kpoints.from_file("KPOINTS")
    #print ("kpoint old",kp_old)
    #print ("kpoint ",kpoints)
    if set(kp_old.kpts[0]) !=set(kpoints.kpts[0]):
                print ("different kp",set(kpoints.kpts[0]),set(kp_old.kpts[0]))
    attempts=0
    #print ("job end time for ", jobname,"=",ticks)
    #print ("job end time for ", jobname,"=",time.localtime(time.time()))
   
    #file='custodian.json'
    #with open(file, 'r') as f:
    #     data = json.load(f)
    #print ("data=",data)
    oszicar = Oszicar("OSZICAR")
    contcar=str(os.getcwd())+str('/')+str('CONTCAR')
    final_str=Structure.from_file(contcar)
    natoms=final_str.composition.num_atoms
    enp=float(oszicar.final_energy)/float(final_str.composition.num_atoms)
    data_cal=[]
    os.chdir('../')
    try:
      data_cal.append({'jobname':jobname,'poscar_initial':mat.as_dict(),'poscar_final':final_str.as_dict(),'incar':incar.as_dict(),'kpoints':kpoints.as_dict(),'final_energy':float(oszicar.final_energy),'contcar':final_str.as_dict()})
    except:
       pass
    json_file=str(jobname)+str('.json')
    f_json=open(json_file,'w')
    f_json.write(json.dumps(data_cal,indent=4,cls=MontyEncoder))
    f_json.close()
    return float(oszicar.final_energy),contcar

    #print (("enp,cont"),enp,final_str)
    #sys.exit()
    ##c.run()
    #for j in all_jobs:
    #        final_energy = j.get_final_energy()
    #        cal_log_new.append({"job": j.as_dict(),
    #                            'job_id': j.job_id,
    #                            "corrections": [],
    #                            'final_energy': final_energy})

    #dumpfn(cal_log_new, chkpt_file, cls=MontyEncoder,
    #           indent=4)
    #f.close()
    #contcar=str(j.vis.name)+str('/')+str('CONTCAR')
    #final_str=Structure.from_file(contcar)
    #natoms=final_str.composition.num_atoms
    #enp=float(final_energy)/float(natoms)
    #return enp,contcar


##########################################################
def check_errors(logfile='vasp.out',timeout=18000):
    errors = []
    #print ("going here 12")
    try:
      run=Vasprun("vasprun.xml")
      if run.converged_electronic == False:
          errors.append("unconverged_electronic")
      if run.converged_ionic == False:
          errors.append("unconverged")
    except:
        pass
    #f=open("OUTCAR",'r')
    try:
       oszicar = Oszicar("OSZICAR")
       ionic_steps=len(oszicar.ionic_steps)
       with open("OUTCAR") as f:
             for line in f:

                 if "NSW" in line:
                     try:
                         #nbands = int(d[-1].strip())
                         nsw= int(str(line.split("=")[1]).split("number of steps")[0])
                         if nsw >1 and ionic_steps==nsw:
                            errors.append("unconverged") 
                     except:
                        pass

                 if "NELM" in line:
                     try:
                         #nbands = int(d[-1].strip())
                         #nsw= int(str(line.split("=")[1]).split("number of steps")[0])
                         nelm= int(str(line.split("=")[1]).split(";")[0])
                         electronic_steps=len(os.electronic_steps[-1])
                         if electronic_steps <nelm:
                             print("Electronically converged")
                         else:
                            errors.append("unconverged_electronic") 
            
                     except:
                        pass

 

    except:
          
         pass
              

    with open(logfile, "r") as f:
        for line in f:
            l = line.strip()

            if "WARNING: Sub-Space-Matrix is not hermitian in" in  line:
                     err="subspacematrix"
                     if err not  in errors:
                        errors.append(err)
            if "Tetrahedron method fails for NKPT<" in  line:
                     err="tet"
                     if err not  in errors:
                        errors.append(err)
            if "Fatal error detecting k-mesh" in  line:
                     err="tet"
                     if err not  in errors:
                        errors.append(err)
            if "Fatal error: unable to match k-point" in  line:
                     err="tet"
                     if err not  in errors:
                        errors.append(err)
            if "Routine TETIRR needs special values" in  line:
                     err="tet"
                     if err not  in errors:
                        errors.append(err)
            if "Tetrahedron method fails for NKPT<" in  line:
                     err="tet"
                     if err not  in errors:
                        errors.append(err)
            if "inverse of rotation matrix was not found (increase" in  line:
                     err="inv_rot_ma"
                     if err not  in errors:
                        errors.append(err)
            if "SYMPREC" in  line:
                     err="inv_rot_ma"
                     if err not  in errors:
                        errors.append(err)
            if "Routine TETIRR needs special values" in  line:
                     err="tetirr"
                     if err not  in errors:
                        errors.append(err)
            if "Could not get correct shift" in  line:
                     err="incorrect_shift"
                     if err not  in errors:
                        errors.append(err)
            if "REAL_OPTLAY: internal error" in  line:
                     err="real_optlay"
                     if err not  in errors:
                        errors.append(err)
            if "REAL_OPT: internal ERROR" in  line:
                     err="real_optlay"
                     if err not  in errors:
                        errors.append(err)
            if "ERROR RSPHER" in  line:
                     err="rspher"
                     if err not  in errors:
                        errors.append(err)
            if "DENTET" in  line:
                     err="dentet"
                     if err not  in errors:
                        errors.append(err)
            if "TOO FEW BAND" in  line:
                     err="too_few_bands"
                     if err not  in errors:
                        errors.append(err)
            if "ERROR: the triple product of the basis vectors" in  line:
                     err="triple_product"
                     if err not  in errors:
                        errors.append(err)
            if "Found some non-integer element in rotation matrix" in  line:
                     err="rot_matrix"
                     if err not  in errors:
                        errors.append(err)
            if "BRIONS problems: POTIM should be increased" in  line:
                     err="brions"
                     if err not  in errors:
                        errors.append(err)
            if "internal error in subroutine PRICEL" in  line:
                     err="pricel"
                     if err not  in errors:
                        errors.append(err)
            if "LAPACK: Routine ZPOTRF failed" in  line:
                     err="zpotrf"
                     if err not  in errors:
                        errors.append(err)
            if "One of the lattice vectors is very long (>50 A), but AMIN" in  line:
                     err="amin"
                     if err not  in errors:
                        errors.append(err)
            if "ZBRENT: fatal internal in" in  line:
                     err="zbrent"
                     if err not  in errors:
                        errors.append(err)
            if "ZBRENT: fatal error in bracketing" in  line:
                     err="zrbent"
                     if err not  in errors:
                        errors.append(err)

            if "ERROR in subspace rotation PSSYEVX" in  line:
                     err="pssyevx"
                     if err not  in errors:
                        errors.append(err)
            if "WARNING in EDDRMM: call to ZHEGV failed" in  line:
                     err="eddrmm"
                     if err not  in errors:
                        errors.append(err)
            if "Error EDDDAV: Call to ZHEGV failed" in  line:
                     err="edddav"
                     if err not  in errors:
                        errors.append(err)
            if "Your FFT grids (NGX,NGY,NGZ) are not sufficient" in  line:
                     err="aliasing_incar"
                     if err not  in errors:
                        errors.append(err)


#    with open(logfile, "r") as f:
#        for line in f:
#            l = line.strip()
#            for err, msgs in VaspErrorHandler.error_msgs.items():
#                for msg in msgs:
#                    if l.find(msg) != -1:
#                        # this checks if we want to run a charged
#                        # computation (e.g., defects) if yes we don't
#                        # want to kill it because there is a change in e-
#                        # density (brmix error)
#                        if err == "brmix" and 'NELECT' in incar:
#                            continue
                        errors.append(err)
    #UNCONVERGED
    st = os.stat('vasp.out')
    if time.time() - st.st_mtime > timeout:
           errors.append('Frozenjob')

    return set(errors)

def check_errorss(logfile='log'):
    errors = []
    with open(logfile, "r") as f:
        for line in f:
            l = line.strip()
            for err, msgs in VaspErrorHandler.error_msgs.items():
                for msg in msgs:
                    if l.find(msg) != -1:
                        # this checks if we want to run a charged
                        # computation (e.g., defects) if yes we don't
                        # want to kill it because there is a change in e-
                        # density (brmix error)
                        if err == "brmix" and 'NELECT' in incar:
                            continue
                        errors.append(err)
    return set(errors)
def check_error(logfile='log'):
    errors = set()
    with open(logfile, "r") as f:
        for line in f:
            l = line.strip()
            for err, msgs in VaspErrorHandler.error_msgs.items():
                for msg in msgs:
                    if l.find(msg) != -1:
                        # this checks if we want to run a charged
                        # computation (e.g., defects) if yes we don't
                        # want to kill it because there is a change in e-
                        # density (brmix error)
                        if err == "brmix" and 'NELECT' in incar:
                            continue
                        errors.add(err)
    return len(errors) > 0


#b1=LA.norm(np.array(mat.lattice.reciprocal_lattice_crystallographic.matrix[0]))
#b2=LA.norm(np.array(mat.lattice.reciprocal_lattice_crystallographic.matrix[1]))
#b3=LA.norm(np.array(mat.lattice.reciprocal_lattice_crystallographic.matrix[2]))
#print b1,b2,b3

#length=20
#n1=max(1,length*b1+0.5)
#n2=max(1,length*b2+0.5)
#n3=max(1,length*b3+0.5)
#print n1,n2,n3
def Auto_Kpoints(mat=None,length=20):
    #kp_file=open("KPOINTS","w")
    #line=str(mat.comment)+'\n'
    #kp_file.write(line)
    #line=str(0)+'\n'
    #kp_file.write(line)
    #line=str("Auto")+'\n'
    #kp_file.write(line)
    #line=str(length)+'\n'
    #kp_file.write(line)
    #kp_file.close()
    #pym_kp=Kpoints.from_file("KPOINTS")

    b1=LA.norm(np.array(mat.structure.lattice.reciprocal_lattice_crystallographic.matrix[0]))
    b2=LA.norm(np.array(mat.structure.lattice.reciprocal_lattice_crystallographic.matrix[1]))
    b3=LA.norm(np.array(mat.structure.lattice.reciprocal_lattice_crystallographic.matrix[2]))

    n1=int(max(1,length*b1+0.5))
    n2=int(max(1,length*b2+0.5))
    n3=int(max(1,length*b3+0.5))
    kpp=Kpoints.gamma_automatic(kpts=(n1,n2,n3))
    return kpp

### kpp=Kpoints.gamma_automatic(kpts=(1,1,1))
### kpp.write_file("KPP")

def converg_encut(encut=500,mat=None):
    en1=-10000
    encut1=encut
    convg_encut1=False
    convg_encut2=False
    
    while  convg_encut2 !=True:
    #while convg_encut1 !=True and  convg_encut2 !=True:
        tol=0.001          #change 0.001
        encut_list=[]
        encut_list.append(encut)
        length=10
        encut1=encut+50
        incar_dict = dict(
            PREC = 'Accurate',
            ENCUT = encut,
            ISMEAR = 0,
            IBRION=2,


            EDIFF = '1E-7',
            NSW = 1,
            NELM = 400,
            NPAR = 4,
            LCHARG = '.FALSE.',
            LWAVE = '.FALSE.' )
        incar = Incar.from_dict(incar_dict)
        kpoints=Auto_Kpoints(mat=mat,length=length)
        print ("running smart_converge for",str(mat.comment)+str('-')+str('ENCUT')+str('-')+str(encut))
        en2,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('ENCUT')+str(mat.comment)+str('-')+str(encut))
        while abs(en2-en1)>tol:
           en1=en2
           encut1=encut+50
           encut_list.append(encut)
           print("Incrementing encut",encut)
           incar_dict = dict(
               PREC = 'Accurate',
               IBRION=2,
               ENCUT = encut1,
               ISMEAR = 0,
               EDIFF = '1E-7',
               NSW = 1,
               NELM = 400,


               NPAR = 4,
               LCHARG = '.FALSE.',
               LWAVE = '.FALSE.' )
           incar = Incar.from_dict(incar_dict)
           print ("running smart_converge for",str(mat.comment)+str('-')+str('ENCUT')+str('-')+str(encut))
           en2,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('ENCUT')+str(mat.comment)+str('-')+str(encut))
        convg_encut1=True


# Some extra points to check
        print ("Some extra points to check for ENCUT")
        
        encut2=encut1+50
        incar['ENCUT']=encut2
        en3,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('ENCUT')+str(mat.comment)+str('-')+str(encut2))
       
        encut3=encut2+50
        incar['ENCUT']=encut3
        en4,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('ENCUT')+str(mat.comment)+str('-')+str(encut3))



        encut4=encut3+50
        incar['ENCUT']=encut4
        en5,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('ENCUT')+str(mat.comment)+str('-')+str(encut4))

        encut5=encut4+50
        incar['ENCUT']=encut5
        en6,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('ENCUT')+str(mat.comment)+str('-')+str(encut5))

        encut6=encut5+50
        incar['ENCUT']=encut6
        en7,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('ENCUT')+str(mat.comment)+str('-')+str(encut6))

        #if en3-en2>tol or en4-en2>tol or en5-en2>tol or en6-en2>tol or en7-en2>tol:
        #if abs(en3-en2)>tol and abs(en4-en2)>tol and abs(en5-en2)>tol and abs(en6-en2)>tol and abs(en7-en2)>tol:
        if abs(en3-en2)>tol or abs(en4-en2)>tol or abs(en5-en2)>tol or abs(en6-en2)>tol or abs(en7-en2)>tol:

           en1=en3
           encut=encut1
           fen=open("EXTRA_ENCUT","w")
           line=str("Extra ENCUT needed ")+str(encut)+'\n'
           fen.write(line)
           fen.close()
        else:
          print ("ENCUT convergence achieved for ",mat.comment,encut)
          convg_encut2=True
    return encut

    
def converg_kpoints(length=0,mat=None):
    en1=-10000
    encut=550
    convg_kp1=False
    convg_kp2=False
    length1=length
    kp_list=[]
    while   convg_kp2 !=True:
    #while convg_kp1 !=True and  convg_kp2 !=True:
        tol=0.001          #change 0.001
        incar_dict = dict(
            PREC = 'Accurate',
            ENCUT = encut,
            ISMEAR = 0,
            IBRION=2,
            EDIFF = '1E-7',
            NSW = 1,
            NELM = 400,

            NPAR = 4,
            LCHARG = '.FALSE.',
            LWAVE = '.FALSE.' )
        incar = Incar.from_dict(incar_dict)
        length1=length1+5
        print ("Incrementing length",length1)
        kpoints=Auto_Kpoints(mat=mat,length=length1)
        mesh=kpoints.kpts[0]
        if mesh not in kp_list:
           kp_list.append(mesh)
           en2,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('KPOINTS')+str(mat.comment)+str('-')+str(length1))



           while abs(en2-en1)>tol:
              en1=en2
              print ("Incrementing length",length1)
              incar_dict = dict(
                  PREC = 'Accurate',
                  IBRION=2,
                  ENCUT = encut,
                  ISMEAR = 0,
                  EDIFF = '1E-7',
                  NSW = 1,

                  NELM = 400,
                  NPAR = 4,
                  LCHARG = '.FALSE.',
                  LWAVE = '.FALSE.' )
              incar = Incar.from_dict(incar_dict)
              while mesh in kp_list:
                 length1=length1+5
                 ##Assuming you are not super unlucky
                 kpoints=Auto_Kpoints(mat=mat,length=length1)
                 mesh=kpoints.kpts[0]
                    
              kpoints=Auto_Kpoints(mat=mat,length=length1)
              mesh=kpoints.kpts[0]
              if mesh not in kp_list:
                 kp_list.append(mesh)
                 en2,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('KPOINTS')+str(mat.comment)+str('-')+str(length1))
              else:
                 length1=length1+5
                 ##Assuming you are not super unlucky
                 kpoints=Auto_Kpoints(mat=mat,length=length1)
                 mesh=kpoints.kpts[0]
                 kp_list.append(mesh)
                 en2,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('KPOINTS')+str(mat.comment)+str('-')+str(length1))
           convg_kp1=True


# Some extra points to check
           print ("Some extra points to check for KPOINTS")
           length3=length1+5
           kpoints=Auto_Kpoints(mat=mat,length=length3)
           mesh=kpoints.kpts[0]
           kp_list.append(mesh)
           en3,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('KPOINTS')+str(mat.comment)+str('-')+str(length3))
       
           length4=length3+5
           kpoints=Auto_Kpoints(mat=mat,length=length4)
           mesh=kpoints.kpts[0]
           kp_list.append(mesh)
           en4,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('KPOINTS')+str(mat.comment)+str('-')+str(length4))


           length5=length4+5
           kpoints=Auto_Kpoints(mat=mat,length=length5)
           mesh=kpoints.kpts[0]
           kp_list.append(mesh)
           en5,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('KPOINTS')+str(mat.comment)+str('-')+str(length5))



           length6=length5+5
           kpoints=Auto_Kpoints(mat=mat,length=length6)
           mesh=kpoints.kpts[0]
           kp_list.append(mesh)
           en6,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('KPOINTS')+str(mat.comment)+str('-')+str(length6))
           length7=length6+5
           kpoints=Auto_Kpoints(mat=mat,length=length7)
           mesh=kpoints.kpts[0]
           kp_list.append(mesh)
           en7,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('KPOINTS')+str(mat.comment)+str('-')+str(length7))



           #if en3-en2>tol or en4-en2>tol or en5-en2>tol or en6-en2>tol or en7-en2>tol:
           #if abs(en3-en2)>tol and abs(en4-en2)>tol and abs(en5-en2)>tol and abs(en6-en2)>tol and abs(en7-en2)>tol:
           if abs(en3-en2)>tol or abs(en4-en2)>tol or abs(en5-en2)>tol or abs(en6-en2)>tol or abs(en7-en2)>tol:
              fkp=open("EXTRA_KPOINTS","w")
              line=str("Extra KPOINTS needed ")+str(length1)+'\n'
              fkp.write(line)
              line=str("en2 length1 ")+str (" ")+str(en2)+str(" ")+str(length1)+'\n'
              fkp.write(line)
              line=str("en3 length3 ")+str (" ")+str(en3)+str(" ")+str(length3)+'\n'
              fkp.write(line)
              line=str("en4 length4 ")+str (" ")+str(en4)+str(" ")+str(length4)+'\n'
              fkp.write(line)
              line=str("en5 length5 ")+str (" ")+str(en5)+str(" ")+str(length5)+'\n'
              fkp.write(line)
              fkp.close()
              en1=en3
              length1=length3
           else:
              print ("KPOINTS convergence achieved for ",mat.comment,length1)
              convg_kp2=True


    return length1

def smart_converge(mat=None,band_str=True,elast_prop=True,optical_prop=True,Raman_calc=False):
    encut=converg_encut(encut=500,mat=mat)
    leng= converg_kpoints(length=0,mat=mat)
    kpoints=Auto_Kpoints(mat=mat,length=leng)
    isif=2
    commen=str(mat.comment)
    lcharg='.FALSE.' 
    if commen.split('@')[0] =='bulk' :
       isif=3
       lcharg='.TRUE.' 
    if commen.split('@')[0] == 'sbulk':
       isif=3
    incar_dict = dict(
        PREC = 'Accurate',
        ENCUT = encut,
        ISMEAR = 0,
        EDIFF = '1E-7',
        EDIFFG = '-1E-3',
        ISIF = 3,


        NEDOS = 5000,
        IBRION = 2,
        NSW = 500,   #change 400
        NELM = 500,  #change 400
        LORBIT = 11,
        LVTOT = '.TRUE.',
        LVHAR = '.TRUE.',
        ISPIN = 2,
        NPAR = 4,
        LCHARG = lcharg,
        LWAVE = '.FALSE.' )
    incar = Incar.from_dict(incar_dict)
    #if commen.startswith('Surf-') :
    #   pol=check_polar(mat_f.structure)
    #   if pol==True:
    #        ase_atoms = AseAtomsAdaptor().get_atoms(mat_f.sttucture)
    #        COM=ase_atoms.get_center_of_mass(scaled=True)
    #        incar.update({"LDIPOL": '.TRUE.',"IDIPOL":4,"ISYM": 0,"DIPOL":COM})
    print ("running smart_converge for",str(mat.comment)+str('-')+str('MAIN-RELAX'))
    en2,contc=run_job(mat=mat,incar=incar,kpoints=kpoints,jobname=str('MAIN-RELAX')+str('-')+str(mat.comment)) 
    cwd=str(os.getcwd()) 
    path=str(contc.split('/CONTCAR')[0])+str('/vasprun.xml')
    v=open(path,"r").readlines()
    for line in v:
           if "NBANDS" in  line:
               nbands=int(line.split(">")[1].split("<")[0])
               print ("nbands=",nbands)
               break
    strt=Structure.from_file(contc)
    mat_f=Poscar(strt)
    mat_f.comment=str(mat.comment)
    if band_str==True:
       incar_dict = dict(
           PREC = 'Accurate',
           ENCUT = encut,
           NBANDS=int(nbands)+10,
           ISMEAR = 0,
           EDIFF = '1E-7',


           LCHARG = '.FALSE.',
           NEDOS = 5000,
           ISPIN = 2,
           ISIF = 2,
           IBRION = 1,
           NELM = 400,
           LORBIT = 11,
           NPAR = 4,
           LWAVE = '.FALSE.' )
       incar = Incar.from_dict(incar_dict)
       user_incar_settings={"EDIFF":1E-6,"ISIF":2,"NSW":0,"LORBIT":11,"ENCUT":encut,"LWAVE":'.FALSE.',"PREC":'Accurate'}
       mpvis = MPNonSCFVaspInputSet(user_incar_settings=incar)


       try: 
           print ("running MAIN-BAND")
           kpoints=mpvis.get_kpoints(mat_f.structure)
           en2B,contcB=run_job(mat=mat_f,incar=incar,kpoints=kpoints,jobname=str('MAIN-BAND')+str('-')+str(mat_f.comment))  
          # kpoints=mpvis.get_kpoints(mat_f.structure)
          # en2B,contcB=run_job(mat=mat_f,incar=incar,kpoints=kpoints,jobname=str('MAIN-BAND')+str('-')+str(mat_f.comment))  
       except:
           print ("No band str calc.")
           if str(os.getcwd)!=cwd:
                print ("Changing directory")
                line=str("cd ")+str(cwd)
                os.chdir(cwd)
                print (os.getcwd())

           pass
    if optical_prop==True:

       incar_dict = dict(
           PREC = 'Accurate',
           ENCUT = encut,
           NBANDS = 3*int(nbands),


           LOPTICS = '.TRUE.',
           ISMEAR = 0,
           EDIFF = '1E-7',

           NEDOS = 5000,
           LCHARG = '.FALSE.',
           ISIF = 2,
           IBRION = 1,
           NELM = 400,
           LORBIT = 11,
           NPAR = 4,
           LWAVE = '.FALSE.' )
       incar = Incar.from_dict(incar_dict)
       kpoints=Auto_Kpoints(mat=mat_f,length=leng)
       try:
           en2OP,contcOP=run_job(mat=mat_f,incar=incar,kpoints=kpoints,jobname=str('MAIN-OPTICS')+str('-')+str(mat_f.comment)) 
       except:
           pass 
    if elast_prop==True:
       incar_dict = dict(
           PREC = 'Accurate',
           ENCUT = float(encut)*1.3,
           ISMEAR = 0,


           ISIF = 3,
           POTIM = 0.015,
           NEDOS = 5000,
           EDIFF = '1E-7',
           IBRION = 6,
           NELM = 400,
           NPAR = 16,
           LWAVE = '.FALSE.' )
       incar = Incar.from_dict(incar_dict)
       sg_mat = SpacegroupAnalyzer(mat_f.structure)
       mat_cvn = sg_mat.get_conventional_standard_structure()
       mat_cvn.sort()
       kpoints=Auto_Kpoints(mat=Poscar(mat_cvn),length=leng)
       try:
          en2E,contcE=run_job(mat=Poscar(mat_cvn),incar=incar,kpoints=kpoints,jobname=str('MAIN-ELASTIC')+str('-')+str(mat_f.comment)) 
       except:
           pass 
    if Raman_calc==True:
        Raman(strt=mat_f,encut=encut,length=leng)
    return en2,mat_f
#####################################################




def bandstr(contc=None,kpoints=None,encut=500):
    #if commen.split('@')[0] =='bulk' :
       incar_dict = dict(
           PREC = 'Accurate',
           ENCUT = encut,

           ISMEAR = 0,
           EDIFF = '1E-4',
           LCHARG = lcharg,
           NEDOS = 5000,
           ISPIN = 2,
           ISIF = 2,
           IBRION = 1,
           NELM = 400,


           LORBIT = 11,
           NPAR = 4,
           LWAVE = '.FALSE.' )
       incar = Incar.from_dict(incar_dict)
       user_incar_settings={"EDIFF":1E-6,"ISIF":2,"NSW":0,"LORBIT":11,"ENCUT":encut,"LWAVE":'.FALSE.',"PREC":'Accurate'}
       mpvis = MPNonSCFVaspInputSet(user_incar_settings=incar)
       final_str=Structure.from_file(contc)
       if kpoints!=None:
          print ("Band structure calculation for user input")
       else:

           kpoints=mpvis.get_kpoints(final_str)
       print ("running smart_converge for",str(mat_f.comment)+str('-')+str('BAND')+str('-')+str(den))
       try:
          en2B,contcB=run_job(mat=mat_f,incar=incar,kpoints=kpoints,jobname=str('MAIN-BAND')+str('-')+str(mat_f.comment))  
       except:
           pass

def elastic_prop(mat_f=None,kpoints=None,encut=500):


       path=str(mat_f.split('/CONTCAR')[0])+str('/vasprun.xml')
       v=open(path,"r").readlines()
       for line in v:
           if "NBANDS" in  line:
               nbands=int(line.split(">")[1].split("<")[0])
               print ("nbands=",nbands)
               break
       incar_dict = dict(
           PREC = 'Accurate',
           ENCUT = encut,
           ISMEAR = 0,
           ISIF = 3,


           POTIM = 0.015,
           NEDOS = 5000,
           EDIFF = '1E-6',
           IBRION = 6,
           NELM = 400,
           NPAR = 16,
           LWAVE = '.FALSE.' )
       incar = Incar.from_dict(incar_dict)
       #v.close()
       final_str=Structure.from_file(contc)
       
       #num_kpoints = (den)*mat.structure.lattice.reciprocal_lattice.volume
       #kpoints = Kpoints.automatic_density(mat.structure, num_kpoints * mat.structure.num_sites,force_gamma=True)
       #kpoints = Kpoints.gamma_automatic(kpts=kpts)
       print ("running smart_converge for",str(mat_f.comment)+str('-')+str('ELASTIC'))
       sg_mat = SpacegroupAnalyzer(mat_f.structure)
       mat_cvn = sg_mat.get_conventional_standard_structure()
       mat_cvn.sort()
       en2E,contcE=run_job(mat=Poscar(mat_cvn),incar=incar,kpoints=kpoints,jobname=str('MAIN-ELASTIC')+str('-')+str(mat_f.comment))  
       #try:
       #   phonts(mat_f)
       #except:
       #   pass





def optical_prop(mat_f=None,kpoints=None,encut=500,nbands=100):
       incar_dict = dict(
           PREC = 'Accurate',
           ENCUT = encut,
           NBANDS = 3*int(nbands),
           LOPTICS = '.TRUE.',
           ISMEAR = 0,
           EDIFF = '1E-4',

           NEDOS = 5000,
           LCHARG = '.FALSE.',
           ISIF = 2,


           IBRION = 1,
           NELM = 400,
           LORBIT = 11,
           NPAR = 4,
           LWAVE = '.FALSE.' )
       incar = Incar.from_dict(incar_dict)
       #num_kpoints = (den)*mat.structure.lattice.reciprocal_lattice.volume
       #kpoints = Kpoints.automatic_density(mat.structure, num_kpoints * mat.structure.num_sites,force_gamma=True)
       kpoints = Kpoints.gamma_automatic(kpts=kpts)
       en2OP,contcOP=run_job(mat=mat_f,incar=incar,kpoints=kpoints,jobname=str('MAIN-OPTICS')+str('-')+str(mat_f.comment))  

def smart_vac(strt=None):
    vac_arr=[]
    sg_mat = SpacegroupAnalyzer(strt)
    mat_cvn = sg_mat.get_conventional_standard_structure()
    mat_cvn.sort()

    cellmax=1 #int(mat_cvn.composition.num_atoms)+int(mat_cvn.ntypesp)#5
    ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
    ase_atoms=ase_atoms*(cellmax,cellmax,cellmax)
    #if len(ase_atoms) >200:
    #   cellmax=1
    #else:
    #   cellmax=2
    #cellmax=int(mat_cvn.composition.num_atoms)+int(mat_cvn.ntypesp)#5
    print ("type of trt is= celmmax",type(strt),cellmax)
    try:
       print ("int(strt.composition.num_atoms)",int(strt.composition.num_atoms))
       print (int(strt.ntypesp))
    except:
       pass
    #cellmax=int(strt.composition.num_atoms)+int(strt.ntypesp)
    vac_done=0
    tol=0.1   #change 0.1
    vac=vac_antisite_def_struct_gen(cellmax=cellmax,struct=strt)
    def_list=[100000  for y in range(len(vac)-1)]
    while vac_done !=1:
        vac=vac_antisite_def_struct_gen(cellmax=cellmax,struct=strt)
        if vac not in vac_arr:
           vac_arr.append(vac)
           print ("in smart_vac(strt=None), cellmax,vac,vac_done=",cellmax,vac,vac_done)
           def_list2,header_list=def_energy(vac=vac)
           diff=matrix(def_list)-matrix(def_list2)
           diff_arr=np.array(diff).flatten()
           print ("in smart_vac(strt=None diff_arr=",diff_arr)
           if any(diff_arr)>tol:
          # for el in diff_arr:
          #     if abs(el)>tol :
                  # print ("in smart_vac(strt=None abs_el=",abs(el))
                   vac_done=0
                   cellmax=cellmax+1
                   ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
                   ase_atoms=ase_atoms*(cellmax,cellmax,cellmax)
                   if len(ase_atoms) >100:
                      vac_done=1
                   def_list=def_list2
           else:
               vac_done=1
#        cellmax=cellmax+1
    return def_list,header_list

def def_energy(vac=[]):
    def_list=[]
    header_list=[]
    fi=str('DEF.INFO')
    f=open(fi,'w')
    for v in vac:
        enp,contc=smart_converge(mat=v)
        #enp,contc=run_job(mat=v,incar=incar,kpoints=kpoints,jobname=str(v.comment))
        strt=Structure.from_file(contc)
        comm=str(v.comment)
        print ("running def_energy for =",comm)
        header_list.append(comm)
        if comm.split('@')[0]=='bulk':
           gs_energy=float(enp)*float(strt.composition.num_atoms)
           print ("in def_energy gs_energy for",comm, gs_energy)
        else:
           chem_pot=sum_chem_pot(strt)
           if comm.startswith("intl"):
              chem_pot=0.0-float(chem_pot)
           def_en=(float(enp)*float(strt.composition.num_atoms))-float(gs_energy)+chem_pot
           print ("in def_energy def_en for",comm, def_en,"chem_pot",chem_pot)
           def_list.append(def_en)
           print (v.comment,'=',def_en)
           line= str(v.comment)+str('=')+str(def_en)+'\n'
           f.write(line)
    f.close()
    return def_list,header_list

def smart_surf(strt=None):
    sg_mat = SpacegroupAnalyzer(strt)
    mat_cvn = sg_mat.get_conventional_standard_structure()
    mat_cvn.sort()
    layers=2
    indices = get_symmetrically_distinct_miller_indices(mat_cvn, 1)
    ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
    for i in indices:
        ase_slab = surface(ase_atoms, i, layers)
        ase_slab.center(vacuum=15, axis=2)
        if len(ase_slab) < 50:
           layers=3
    surf_arr=[]
    surf_done=0
    tol=0.1 #change 0.01
    surf=surfer(mat=strt,layers=layers)
    surf_list=[100000  for y in range(len(surf)-1)]
    print ("in smart_surf :surf,surf_list=",surf, surf_list)
    while surf_done !=1:
        layers=layers+1
        indices = get_symmetrically_distinct_miller_indices(mat_cvn, 1)
        ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
        for i in indices:
            ase_slab = surface(ase_atoms, i, layers)
            ase_slab.center(vacuum=15, axis=2)
            if len(ase_slab) > 100:
               surf_done=1
            if (ase_slab.get_cell()[2][2]) > 40:
               surf_done=1
        surf=surfer(mat=strt,layers=layers)
        if surf not in surf_arr:
           surf_arr.append(surf)
           surf_list2,surf_header_list=surf_energy(surf=surf)
           print ("in smart_surf :surf2,surf_list2=",surf_list2,surf_header_list)
           diff=matrix(surf_list)-matrix(surf_list2)
           print ("in smart_surf :surf3,surf_list3=",matrix(surf_list),matrix(surf_list2))
           diff_arr=np.array(diff).flatten()
           if any(diff_arr)>tol:
           #for el in diff_arr:
           #    if abs(el)>tol :
           #        print ("in smart_surf :abs el=",abs(el))
                   surf_done=0
                   surf_list=surf_list2
           else:
                surf_done=1
    return surf_list,surf_header_list
def surf_energy(surf=[]):
    surf_list=[]
    surf_header_list=[]
    fi=str('SURF.INFO')
    f=open(fi,'w')
    for s in surf:
        enp,contc=smart_converge(mat=s)
        strt=Structure.from_file(contc)
        m=strt.lattice.matrix
        if s.comment.split('@')[0]=='sbulk':

           gs_energy=enp
           print ("in surf_energy gs_energy for",s.comment, gs_energy)
        else:
           surf_en=16*(float(enp)*(strt.composition.num_atoms)-float(strt.composition.num_atoms)*float(gs_energy))/(2*np.linalg.norm(np.cross(m[0], m[1])))
           print ("in surf_energy surf_en for",s.comment, surf_en)
           surf_list.append(surf_en)
           print (s.comment,'=',surf_en)
           line= str(s.comment)+str('=')+str(surf_en)+'\n'
           f.write(line)
    f.close()
    return surf_list,surf_header_list

def make_big(poscar=None):
    struct=poscar.structure
    comm=poscar.comment
    a, b, c = struct.lattice.abc
    size=11
    struct.make_supercell([int(float(size)/float(a))+1,int(float(size)/float(b))+1,int(float(size)/float(c))+1])
    big=Poscar(struct)
    big.comment=str(comm)
    return big


def get_smart_surf_def(mat=None):
    #with MPRester() as mp:
    #     strt = mp.get_structure_by_material_id(mpid)
    mat=Poscar.from_file("POSCAR")
#
    #     mpid=mpid.replace('-','_')
    #     mpid=str('bulk@')+str(mpid)
    #     mat.comment=mpid

    print ("Running Main Relax")
    en,contc=smart_converge(mat=mat)
    relaxed_str=Poscar.from_file(contc)
    commen=str(relaxed_str.comment)
    kpfile=str("MAIN-RELAX-")+str(commen)+str('/KPOINTS')
    kp=Kpoints.from_file(kpfile)
    incfile=str("MAIN-RELAX-")+str(commen)+str('/INCAR')
    inc=Incar.from_file(incfile)
    encut=inc['ENCUT']

    #pos1=Poscar.from_file("POSCAR")
    #kp1=Kpoints.from_file("KPOINTS")
    #print kp1.kpts
    tag=True
    kppa=1
    while tag:
      kpoints = Kpoints.automatic_density(relaxed_str.structure, kppa,force_gamma=True)
      #print kpoints.kpts[0],kp1.kpts[0]
      kppa=kppa+1
      if (kpoints.kpts[0])>=(kp1.kpts[0]):
         print (kpoints.kpts,kp1.kpts)
         tag=False
      print ("kpp=",kpp)

    kpoints = Kpoints.automatic_density(relaxed_str.structure, kppa,force_gamma=True)
    elastic_prop(mat_f=relaxed_str,kpoints=kpoints,encut=encut)


    bandstr(contc=relaxed_str,encut=encut)
    optical_prop(mat_f=relaxed_str,kpoints=kpoints,encut=encut,nbands=40)
    big=make_big(relaxed_str)
    try:
        Raman(strt=big,encut=encut,kpoints=kp)
    except:
        print ("Failed Raman")
        pass
#get_smart_surf_def(mat='POSCAR')
def main_func(mpid='',mat=None):
    if mpid !='':
       with MPRester() as mp:
         strt = mp.get_structure_by_material_id(mpid)
         sg_mat = SpacegroupAnalyzer(strt)
         mat_cvn = sg_mat.get_conventional_standard_structure()
         mat_cvn.sort()
         if int(strt.composition._natoms)==int(mat_cvn.composition._natoms):
             mat= Poscar(mat_cvn)
         else:
             mat=Poscar(strt)



         mpid=mpid.replace('-','_')
         mpid=str('bulk@')+str(mpid)
         mat.comment=mpid

    en,final=smart_converge(mat=mat)
    print (en,contc)
#main_func(mpid='mp-782')
