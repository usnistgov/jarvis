"""
Module to run OptB88vdW based High-throughput calculations
"""
from __future__ import division, unicode_literals, print_function
import os, socket, shutil
from custodian.vasp.jobs import VaspJob
from pymatgen.io.vasp import VaspInput, Vasprun
from pymatgen.io.vasp.outputs import Oszicar
from subprocess import Popen, PIPE
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import sys, shutil, glob, codecs
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import (
    Slab,
    SlabGenerator,
    generate_all_slabs,
    get_symmetrically_distinct_miller_indices,
)
from custodian.vasp.handlers import (
    VaspErrorHandler,
    UnconvergedErrorHandler,
    MeshSymmetryErrorHandler,
    NonConvergingErrorHandler,
    PotimErrorHandler,
)
from pymatgen.symmetry.bandstructure import HighSymmKpath
import json, yaml
from numpy import linalg as LA
import time
from collections import OrderedDict
from jarvis.lammps.jlammps import vac_antisite_def_struct_gen, surfer

# from jarvis.lammps.Surf_Def import vac_antisite_def_struct_gen,surfer
from pymatgen.ext.matproj import MPRester
import subprocess
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar, Poscar, VaspInput, Potcar, Kpoints

# from pymatgen.io.vasp.sets import MPVaspInputSet, MPNonSCFVaspInputSet
from numpy import matrix
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
import operator
from monty.serialization import loadfn, dumpfn
from monty.json import MontyEncoder, MontyDecoder

ncores = 16
try:
    # Define in enenv_variable file
    main_exe = os.environ["vasp_bulk_exe"]
    surf_exe = os.environ["vasp_surf_exe"]
    nw_exe = os.environ["vasp_nw_exe"]
    soc_exe = os.environ["vasp_soc_exe"]
    pot_yaml = str(os.path.join(os.path.dirname(__file__), "Special_POTCAR.yaml"))
    vdw_dat = os.environ["vasp_vdw_dat"]
    json_dat = os.environ["data_json"]
    ncores = int(os.environ["ncores"])
    nnodes = int(os.environ["nnodes"])
    mem = str(os.environ["mem"])
    walltime = str(os.environ["walltime"])
    mp_cmd = str(os.environ["mp_cmd"])
except:
    pass

optb88dict = dict(
    PREC="Accurate",
    ISMEAR=0,
    IBRION=2,
    GGA="BO",
    PARAM1=0.1833333333,
    PARAM2=0.2200000000,
    LUSE_VDW=".TRUE.",
    AGGAC=0.0000,
    EDIFF="1E-7",
    NSW=1,
    NELM=400,
    ISIF=2,
    NPAR=ncores,
    #NPAR=np.sqrt(ncores),
    LCHARG=".FALSE.",
    LWAVE=".FALSE.",
)
functional = "PBE"
use_incar_dict = optb88dict

# Change use_incar_dict based on a functional,
# examples given for LDA and PBE

# functional='LDA'
# ldadict = dict(
#            PREC = 'Accurate',
#            ISMEAR = 0,
#            IBRION=2,
#
#            EDIFF = '1E-7',
#            NSW = 1,
#            NELM = 400,
#            ISIF = 2,
#            NPAR = np.sqrt(ncores),
#            LCHARG = '.FALSE.',
#            LWAVE = '.FALSE.' )
#
# use_incar_dict=ldadict
#
# functional='PBE'
# pbedict = dict(
#            PREC = 'Accurate',
#            ISMEAR = 0,
#            IBRION=2,
#
#            GGA = 'PE',
#
#            EDIFF = '1E-7',
#            NSW = 1,
#            NELM = 400,
#            ISIF = 2,
#            NPAR = np.sqrt(ncores),
#            LCHARG = '.FALSE.',
#            LWAVE = '.FALSE.' )
#
# use_incar_dict=pbedict


def check_polar(file):
    """
    Check if the surface structure is polar
    by comparing atom types at top and bottom

    Args:
         file:Structure object (surface with vacuum)
    Returns:
           polar:True/False   
    """
    up = 0
    dn = 0
    coords = np.array(file.frac_coords)
    z_max = max(coords[:, 2])
    z_min = min(coords[:, 2])
    for site in file:
        if site.frac_coords[2] == z_max:
            up = up + site.specie.number
        if site.frac_coords[2] == z_min:
            dn = dn + site.specie.number
    polar = False
    if up != dn:
        print("polar")
        polar = True
    if up == dn:
        print("Non-polar")
        polar = False
    return polar


def main_outcar(fil=""):
    cnvg = False
    try:
        f = open(fil, "r")
        lines = f.read().splitlines()
        f.close()
        cnvg = False
        for i in lines:
            if "General timing and accounting informations for this job" in i:
                cnvg = True
        # print fil,cnvg
    except:
        pass
    return cnvg


def get_lowest_en_from_mp(formula, MAPI_KEY="", all_structs=False):
    """
    Lowest energy/chemical potential of an element
    from the materialsproject/jarvis database.
    Note: Get the api key from materialsproject/jarvis website.
 
    Args:
        formula: say Al, Ni etc.
        MAPI_KEY: should be defines in the environment
        all_structs: all structures or just stable ones (True/False)
    Returns:
         enp: energy per atom

    """
    if not MAPI_KEY:
        MAPI_KEY = os.environ.get("MAPI_KEY", "")
        if not MAPI_KEY:
            print("API key not provided")
            print(
                "get API KEY from materialsproject and set it to the MAPI_KEY environment variable. aborting ... "
            )
            sys.exit()
    with MPRester(MAPI_KEY) as m:
        data = m.get_data(formula)
        structures = []
        x = {}
        print(
            "\nnumber of structures matching the chemical formula {0} = {1}".format(
                formula, len(data)
            )
        )
        print(
            "The one with the the lowest energy above the hull is returned, unless all_structs is set to True"
        )
        for d in data:
            mpid = str(d["material_id"])
            x[mpid] = d["e_above_hull"]
            if all_structs:
                structure = m.get_structure_by_material_id(mpid)
                structures.append(structure)
        else:
            mineah_key = sorted(x.items(), key=operator.itemgetter(1))[0][0]
            print(
                "The id of the material corresponding to the lowest energy above the hull = {0}".format(
                    mineah_key
                )
            )
            if mineah_key:
                with MPRester(MAPI_KEY) as m:
                    data = m.get_data(mineah_key)
                    x = {}
                    for d in data:
                        x["energy_per_atom"] = str(d["energy_per_atom"])
                        enp = x["energy_per_atom"]
                # return m.get_structure_by_material_id(mineah_key)
                return enp
            else:
                return None


def sum_chem_pot(strt=None):
    """
    Helper function for sump of chemical potential

    Args:
        strt: Structure object
    Returns:
           sum: sum of energy
    """
    sum = 0
    symb = strt.symbol_set
    for el in symb:
        enp = get_lowest_en_from_mp(el)
        sum = float(sum) + float(enp)
    return sum


def run_job(mat=None, incar=None, kpoints=None, jobname="", copy_file=[]):
    """
    Generic function to run a VASP job, error correction implemented using
    custodian package
    A jobname+.json file is produced after successful completion of the job
    A first_cust.py file is generated which is invoked using python command

    Args:
        mat: Poscar object with structure information
        incar: Incar object with control information
        kpoints: Kpoints object with mesh information
        jobname: a uniq name for a job-type, say MAIN-ELAST for elastic properties
        copy_file: copy file from previous runs, say CHGCAR for Non-SCF bandstructure calculations 
    Returns:
        f_energy: final energy
        contcar: path to final relaxed structure
    """
    # hostname=str(socket.gethostname())
    poscar_list = [(mat)]
    cwd = str(os.getcwd())
    job_name = str(mat.comment)
    job_dir = str(jobname)
    run_file = str(os.getcwd()) + str("/") + str(jobname) + str(".json")
    run_dir = str(os.getcwd()) + str("/") + str(jobname)
    if mat.comment.startswith("Surf"):
        [a, b, c] = kpoints.kpts[0]
        kpoints.kpts = [[a, b, 1]]
        try:
            pol = check_polar(mat.structure)
            if pol == True:
                ase_atoms = AseAtomsAdaptor().get_atoms(mat.structure)
                COM = ase_atoms.get_center_of_mass(scaled=True)
                print("COM=", COM)
                print("Found polar surface, will be setting dipole corrections")
                incar.update(
                    {
                        "LDIPOL": ".TRUE.",
                        "IDIPOL": 3,
                        "ISYM": 0,
                        "DIPOL": str(COM[0])
                        + str(" ")
                        + str(COM[2])
                        + str(" ")
                        + str(COM[2]),
                    }
                )
                print("Polar surface encountered in run_job", mat.comment)
        except:
            pass
    wait = False
    json_file = str(jobname) + str(".json")
    print("json should be here=", str(os.getcwd()) + str("/") + str(json_file))
    # print ('json should be=',json_file,run_file,os.getcwd())
    if os.path.exists(str(os.getcwd()) + str("/") + str(json_file)):
        try:
            data_cal = loadfn(
                str(os.getcwd()) + str("/") + str(json_file), cls=MontyDecoder
            )
            tmp_outcar = (
                str(os.getcwd())
                + str("/")
                + str(json_file.split(".json")[0])
                + str("/OUTCAR")
            )
            print("outcar is", tmp_outcar)
            wait = main_outcar(tmp_outcar)  # True
            print("outcar status", wait)
            if wait == True:
                f_energy = data_cal[0]["final_energy"]
                contcar = (
                    str(os.getcwd())
                    + str("/")
                    + str(json_file.split(".json")[0])
                    + str("/CONTCAR")
                )
                return f_energy, contcar
        except:
            pass
    while wait == False:
        print("Setting up POTCAR")
        with open(pot_yaml, "r") as f:
            doc = yaml.load(f)
            pots = doc["POTCAR"]
        new_symb = []
        for el in mat.site_symbols:
            new_symb.append(pots[el])
        # potcar = Potcar(symbols=new_symb,functional=functional)
        try:
            potcar = Potcar(symbols=new_symb, functional=functional)
        except:
            print(
                "JARVIS-ERROR: Could not set POTCAR, check POTCAR yaml file and VASP_PSP_DIR"
            )
            pass
        if not os.path.exists(run_dir):
            print("Starting new job")
            os.makedirs(run_dir)
            os.chdir(run_dir)
            incar.write_file("INCAR")
            potcar.write_file("POTCAR")
            kpoints.write_file("KPOINTS")
            mat.write_file("POSCAR")
            for i in copy_file:
                print("copying", i)
                shutil.copy2(i, "./")
            # f=open('job.out','w')
            cmd = str("python  first_cust.py >out_dat") + "\n"
            cust_file = open("first_cust.py", "w")
            cline = (
                str(
                    "from pymatgen.io.vasp.inputs import Incar, Poscar, VaspInput,Potcar, Kpoints"
                )
                + "\n"
            )
            cust_file.write(cline)
            cline = str("import os,shutil") + "\n"
            cust_file.write(cline)
            cline = str("from custodian.vasp.jobs import VaspJob") + "\n"
            cust_file.write(cline)
            cline = (
                str(
                    "from custodian.vasp.handlers import VaspErrorHandler, UnconvergedErrorHandler,MeshSymmetryErrorHandler, NonConvergingErrorHandler, PotimErrorHandler"
                )
                + "\n"
            )
            cust_file.write(cline)
            cline = (
                str(
                    "from custodian.vasp.validators import VaspFilesValidator,VasprunXMLValidator"
                )
                + "\n"
            )
            cust_file.write(cline)
            cline = str("from custodian.custodian import Custodian") + "\n"
            cust_file.write(cline)
            cline = str('inc=Incar.from_file("INCAR")') + "\n"
            cust_file.write(cline)
            cline = str('pot=Potcar.from_file("POTCAR")') + "\n"
            cust_file.write(cline)
            cline = str('pos=Poscar.from_file("POSCAR")') + "\n"
            cust_file.write(cline)
            cline = str('kp=Kpoints.from_file("KPOINTS")') + "\n"
            cust_file.write(cline)
            cline = str("shutil.copy2('") + vdw_dat + str("','./')") + "\n"
            cust_file.write(cline)
            cline = str('vinput = VaspInput.from_directory(".")') + "\n"
            cust_file.write(cline)
            cline = (
                str("job=VaspJob(['mpirun', '")
                + str(main_exe)
                + str("'], final=False, backup=False)")
                + "\n"
            )
            if mat.comment.startswith("Surf"):
                cline = (
                    str("job=VaspJob(['mpirun',  '")
                    + str(surf_exe)
                    + str("'], final=False, backup=False)")
                    + "\n"
                )
            if "SOC" in jobname:
                cline = (
                    str("job=VaspJob(['mpirun',  '")
                    + str(soc_exe)
                    + str("'], final=False, backup=False)")
                    + "\n"
                )
            cust_file.write(cline)
            cline = (
                str(
                    "handlers = [VaspErrorHandler(), MeshSymmetryErrorHandler(),UnconvergedErrorHandler(), NonConvergingErrorHandler(),PotimErrorHandler()]"
                )
                + "\n"
            )
            cust_file.write(cline)
            cline = str("validators = [VasprunXMLValidator()]") + "\n"
            # cline=str('validators = [VaspFilesValidator()]')+'\n'
            cust_file.write(cline)
            cline = (
                str("c = Custodian(handlers, [job],max_errors=5,validators=validators)")
                + "\n"
            )
            cust_file.write(cline)
            cline = str("c.run()") + "\n"
            cust_file.write(cline)
            cust_file.close()
            print("I AM HERE 2")
            os.system(cmd)
            if os.path.isfile("OUTCAR"):
                try:
                    wait = main_outcar("OUTCAR")  # Vasprun("vasprun.xml").converged
                except:
                    pass
            print("End of the first loop", os.getcwd(), wait)

        else:
            print("Jobs seens to have started before")
            os.chdir(run_dir)
            wait = False
            if os.path.isfile("OUTCAR"):
                try:
                    wait = main_outcar("OUTCAR")  # Vasprun("vasprun.xml").converged
                    # wait=Vasprun("vasprun.xml").converged
                except:
                    pass
            print("Tried to find OUTCAR, wait is=", wait)
            if wait == False:
                incar.write_file("INCAR")
                kpoints.write_file("KPOINTS")
                mat.write_file("POSCAR")
                try:
                    potcar.write_file("POTCAR")
                    print("FOUND OLD CONTCAR in", os.getcwd())
                    old_contcar = Poscar.from_file("CONTCAR")
                    # old_contcar.write_file('POSCAR')
                    copy_cmd = str("cp CONTCAR POSCAR")
                    print("copy_cmd=", copy_cmd)
                    if "ELAST" not in jobname:
                        # Because in ELASTIC calculations structures are deformed
                        os.system(copy_cmd)
                    # time.sleep(3)
                except:
                    pass
                for i in copy_file:
                    print("copying", i)
                    shutil.copy2(i, "./")

                cmd = str("python  first_cust.py >out_dat") + "\n"
                cust_file = open("first_cust.py", "w")
                cline = (
                    str(
                        "from pymatgen.io.vasp.inputs import Incar, Poscar, VaspInput,Potcar, Kpoints"
                    )
                    + "\n"
                )
                cust_file.write(cline)
                cline = str("import os,shutil") + "\n"
                cust_file.write(cline)
                cline = str("from custodian.vasp.jobs import VaspJob") + "\n"
                cust_file.write(cline)
                cline = (
                    str(
                        "from custodian.vasp.handlers import VaspErrorHandler, UnconvergedErrorHandler,MeshSymmetryErrorHandler, NonConvergingErrorHandler, PotimErrorHandler"
                    )
                    + "\n"
                )
                cust_file.write(cline)
                cline = (
                    str(
                        "from custodian.vasp.validators import VaspFilesValidator,VasprunXMLValidator"
                    )
                    + "\n"
                )
                cust_file.write(cline)
                cline = str("from custodian.custodian import Custodian") + "\n"
                cust_file.write(cline)
                cline = str('inc=Incar.from_file("INCAR")') + "\n"
                cust_file.write(cline)
                cline = str('pot=Potcar.from_file("POTCAR")') + "\n"
                cust_file.write(cline)
                cline = str('pos=Poscar.from_file("POSCAR")') + "\n"
                cust_file.write(cline)
                cline = str('kp=Kpoints.from_file("KPOINTS")') + "\n"
                cust_file.write(cline)
                cline = str("shutil.copy2('") + vdw_dat + str("','./')") + "\n"
                cust_file.write(cline)
                cline = str('vinput = VaspInput.from_directory(".")') + "\n"
                cust_file.write(cline)
                cline = (
                    str("job=VaspJob(['mpirun', '")
                    + str(main_exe)
                    + str("'], final=False, backup=False)")
                    + "\n"
                )
                if mat.comment.startswith("Surf"):
                    cline = (
                        str("job=VaspJob(['mpirun',  '")
                        + str(surf_exe)
                        + str("'], final=False, backup=False)")
                        + "\n"
                    )
                if "SOC" in jobname:
                    cline = (
                        str("job=VaspJob(['mpirun',  '")
                        + str(soc_exe)
                        + str("'], final=False, backup=False)")
                        + "\n"
                    )
                cust_file.write(cline)
                cline = (
                    str(
                        "handlers = [VaspErrorHandler(), MeshSymmetryErrorHandler(),UnconvergedErrorHandler(), NonConvergingErrorHandler(),PotimErrorHandler()]"
                    )
                    + "\n"
                )
                cust_file.write(cline)
                # cline=str('validators = [VaspFilesValidator()]')+'\n'
                cline = str("validators = [VasprunXMLValidator()]") + "\n"
                cust_file.write(cline)
                cline = (
                    str(
                        "c = Custodian(handlers, [job],max_errors=5,validators=validators)"
                    )
                    + "\n"
                )
                cust_file.write(cline)
                cline = str("c.run()") + "\n"
                cust_file.write(cline)
                cust_file.close()
                os.system(cmd)
                if os.path.isfile("OUTCAR"):
                    try:
                        wait = main_outcar("OUTCAR")  # Vasprun("vasprun.xml").converged
                        # wait=Vasprun("vasprun.xml").converged
                    except:
                        pass
    f_energy = "na"
    enp = "na"
    contcar = str(os.getcwd()) + str("/") + str("CONTCAR")
    final_str = Structure.from_file(contcar)
    try:
        oszicar = Oszicar("OSZICAR")
        f_energy = float(oszicar.final_energy)
        enp = float(oszicar.final_energy) / float(final_str.composition.num_atoms)
    except:
        print("Error in OSZICAR file during re-run jpb")
        pass
    natoms = final_str.composition.num_atoms
    os.chdir("../")
    if wait == True:
        data_cal = []
        data_cal.append(
            {
                "jobname": jobname,
                "poscar_initial": mat.as_dict(),
                "poscar_final": final_str.as_dict(),
                "incar": incar.as_dict(),
                "kpoints": kpoints.as_dict(),
                "final_energy": (f_energy),
                "contcar": final_str.as_dict(),
            }
        )
        json_file = str(jobname) + str(".json")
        f_json = open(json_file, "w")
        f_json.write(json.dumps(data_cal, indent=4, cls=MontyEncoder))
        f_json.close()
        print("Wrote json file", contcar)
        return f_energy, contcar


def Auto_Kpoints(mat=None, length=20):
    """
    Geting Kpoints object from structure and line-density

    Args:
         mat: Poscar object with structure information
         length: line-density
    Returns:
         kpp: Kpoint object
    """

    b1 = LA.norm(
        np.array(mat.structure.lattice.reciprocal_lattice_crystallographic.matrix[0])
    )
    b2 = LA.norm(
        np.array(mat.structure.lattice.reciprocal_lattice_crystallographic.matrix[1])
    )
    b3 = LA.norm(
        np.array(mat.structure.lattice.reciprocal_lattice_crystallographic.matrix[2])
    )

    n1 = int(max(1, length * b1 + 0.5))
    n2 = int(max(1, length * b2 + 0.5))
    n3 = int(max(1, length * b3 + 0.5))
    kpp = Kpoints.gamma_automatic(kpts=(n1, n2, n3))
    return kpp


def converg_encut(encut=500, mat=None):
    """
    Function to converg plane-wave cut-off

    Args:
        encut: intial cutoff
        mat: Poscar object
    Returns:
           encut: converged cut-off
    """
    en1 = -10000
    encut1 = encut
    convg_encut1 = False
    convg_encut2 = False

    while convg_encut2 != True:
        # while convg_encut1 !=True and  convg_encut2 !=True:
        tol = 0.001  # change 0.001
        encut_list = []
        encut_list.append(encut)
        length = 10
        encut1 = encut + 50
        incar_dict = use_incar_dict
        incar_dict.update({"ENCUT": encut})
        # print (use_incar_dict)
        incar = Incar.from_dict(incar_dict)
        kpoints = Auto_Kpoints(mat=mat, length=length)
        print(
            "running smart_converge for",
            str(mat.comment) + str("-") + str("ENCUT") + str("-") + str(encut),
        )
        en2, contc = run_job(
            mat=mat,
            incar=incar,
            kpoints=kpoints,
            jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut),
        )
        while abs(en2 - en1) > tol:
            en1 = en2
            encut1 = encut + 50
            encut_list.append(encut)
            print("Incrementing encut", encut)
            incar_dict = use_incar_dict
            incar_dict.update({"ENCUT": encut1})
            incar = Incar.from_dict(incar_dict)
            print(
                "running smart_converge for",
                str(mat.comment) + str("-") + str("ENCUT") + str("-") + str(encut),
            )
            en2, contc = run_job(
                mat=mat,
                incar=incar,
                kpoints=kpoints,
                jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut),
            )
        convg_encut1 = True

        # Some extra points to check
        print("Some extra points to check for ENCUT")

        encut2 = encut1 + 50
        incar["ENCUT"] = encut2
        en3, contc = run_job(
            mat=mat,
            incar=incar,
            kpoints=kpoints,
            jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut2),
        )

        encut3 = encut2 + 50
        incar["ENCUT"] = encut3
        en4, contc = run_job(
            mat=mat,
            incar=incar,
            kpoints=kpoints,
            jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut3),
        )

        encut4 = encut3 + 50
        incar["ENCUT"] = encut4
        en5, contc = run_job(
            mat=mat,
            incar=incar,
            kpoints=kpoints,
            jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut4),
        )

        encut5 = encut4 + 50
        incar["ENCUT"] = encut5
        en6, contc = run_job(
            mat=mat,
            incar=incar,
            kpoints=kpoints,
            jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut5),
        )

        encut6 = encut5 + 50
        incar["ENCUT"] = encut6
        en7, contc = run_job(
            mat=mat,
            incar=incar,
            kpoints=kpoints,
            jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut6),
        )

        # if en3-en2>tol or en4-en2>tol or en5-en2>tol or en6-en2>tol or en7-en2>tol:
        # if abs(en3-en2)>tol and abs(en4-en2)>tol and abs(en5-en2)>tol and abs(en6-en2)>tol and abs(en7-en2)>tol:
        if (
            abs(en3 - en2) > tol
            or abs(en4 - en2) > tol
            or abs(en5 - en2) > tol
            or abs(en6 - en2) > tol
            or abs(en7 - en2) > tol
        ):

            en1 = en3
            encut = encut1
            fen = open("EXTRA_ENCUT", "w")
            line = str("Extra ENCUT needed ") + str(encut) + "\n"
            fen.write(line)
            fen.close()
        else:
            print("ENCUT convergence achieved for ", mat.comment, encut)
            convg_encut2 = True
    return encut


def converg_kpoints(length=0, mat=None):
    """
    Function to converg K-points

    Args:
        lenght: K-point line density
        mat: Poscar object with structure information
    Returns:
           length1: K-point line density
    """

    en1 = -10000
    encut = 550
    convg_kp1 = False
    convg_kp2 = False
    length1 = length
    kp_list = []
    while convg_kp2 != True:
        # while convg_kp1 !=True and  convg_kp2 !=True:
        tol = 0.001  # change 0.001
        incar_dict = use_incar_dict
        incar_dict.update({"ENCUT": encut})
        incar = Incar.from_dict(incar_dict)
        length1 = length1 + 5
        print("Incrementing length", length1)
        kpoints = Auto_Kpoints(mat=mat, length=length1)
        mesh = kpoints.kpts[0]
        if mesh not in kp_list:
            kp_list.append(mesh)
            en2, contc = run_job(
                mat=mat,
                incar=incar,
                kpoints=kpoints,
                jobname=str("KPOINTS") + str(mat.comment) + str("-") + str(length1),
            )

            while abs(en2 - en1) > tol:
                en1 = en2
                print("Incrementing length", length1)
                incar_dict = use_incar_dict
                incar_dict.update({"ENCUT": encut})
                incar = Incar.from_dict(incar_dict)
                while mesh in kp_list:
                    length1 = length1 + 5
                    ##Assuming you are not super unlucky
                    kpoints = Auto_Kpoints(mat=mat, length=length1)
                    mesh = kpoints.kpts[0]

                kpoints = Auto_Kpoints(mat=mat, length=length1)
                mesh = kpoints.kpts[0]
                if mesh not in kp_list:
                    kp_list.append(mesh)
                    en2, contc = run_job(
                        mat=mat,
                        incar=incar,
                        kpoints=kpoints,
                        jobname=str("KPOINTS")
                        + str(mat.comment)
                        + str("-")
                        + str(length1),
                    )
                else:
                    length1 = length1 + 5
                    ##Assuming you are not super unlucky
                    kpoints = Auto_Kpoints(mat=mat, length=length1)
                    mesh = kpoints.kpts[0]
                    kp_list.append(mesh)
                    en2, contc = run_job(
                        mat=mat,
                        incar=incar,
                        kpoints=kpoints,
                        jobname=str("KPOINTS")
                        + str(mat.comment)
                        + str("-")
                        + str(length1),
                    )
            convg_kp1 = True

            # Some extra points to check
            print("Some extra points to check for KPOINTS")
            length3 = length1 + 5
            kpoints = Auto_Kpoints(mat=mat, length=length3)
            mesh = kpoints.kpts[0]
            kp_list.append(mesh)
            en3, contc = run_job(
                mat=mat,
                incar=incar,
                kpoints=kpoints,
                jobname=str("KPOINTS") + str(mat.comment) + str("-") + str(length3),
            )

            length4 = length3 + 5
            kpoints = Auto_Kpoints(mat=mat, length=length4)
            mesh = kpoints.kpts[0]
            kp_list.append(mesh)
            en4, contc = run_job(
                mat=mat,
                incar=incar,
                kpoints=kpoints,
                jobname=str("KPOINTS") + str(mat.comment) + str("-") + str(length4),
            )

            length5 = length4 + 5
            kpoints = Auto_Kpoints(mat=mat, length=length5)
            mesh = kpoints.kpts[0]
            kp_list.append(mesh)
            en5, contc = run_job(
                mat=mat,
                incar=incar,
                kpoints=kpoints,
                jobname=str("KPOINTS") + str(mat.comment) + str("-") + str(length5),
            )

            length6 = length5 + 5
            kpoints = Auto_Kpoints(mat=mat, length=length6)
            mesh = kpoints.kpts[0]
            kp_list.append(mesh)
            en6, contc = run_job(
                mat=mat,
                incar=incar,
                kpoints=kpoints,
                jobname=str("KPOINTS") + str(mat.comment) + str("-") + str(length6),
            )
            length7 = length6 + 5
            kpoints = Auto_Kpoints(mat=mat, length=length7)
            mesh = kpoints.kpts[0]
            kp_list.append(mesh)
            en7, contc = run_job(
                mat=mat,
                incar=incar,
                kpoints=kpoints,
                jobname=str("KPOINTS") + str(mat.comment) + str("-") + str(length7),
            )

            # if en3-en2>tol or en4-en2>tol or en5-en2>tol or en6-en2>tol or en7-en2>tol:
            # if abs(en3-en2)>tol and abs(en4-en2)>tol and abs(en5-en2)>tol and abs(en6-en2)>tol and abs(en7-en2)>tol:
            if (
                abs(en3 - en2) > tol
                or abs(en4 - en2) > tol
                or abs(en5 - en2) > tol
                or abs(en6 - en2) > tol
                or abs(en7 - en2) > tol
            ):
                fkp = open("EXTRA_KPOINTS", "w")
                line = str("Extra KPOINTS needed ") + str(length1) + "\n"
                fkp.write(line)
                line = (
                    str("en2 length1 ")
                    + str(" ")
                    + str(en2)
                    + str(" ")
                    + str(length1)
                    + "\n"
                )
                fkp.write(line)
                line = (
                    str("en3 length3 ")
                    + str(" ")
                    + str(en3)
                    + str(" ")
                    + str(length3)
                    + "\n"
                )
                fkp.write(line)
                line = (
                    str("en4 length4 ")
                    + str(" ")
                    + str(en4)
                    + str(" ")
                    + str(length4)
                    + "\n"
                )
                fkp.write(line)
                line = (
                    str("en5 length5 ")
                    + str(" ")
                    + str(en5)
                    + str(" ")
                    + str(length5)
                    + "\n"
                )
                fkp.write(line)
                fkp.close()
                en1 = en3
                length1 = length3
            else:
                print("KPOINTS convergence achieved for ", mat.comment, length1)
                convg_kp2 = True

    return length1


def smart_converge(
    mat=None,
    encut="",
    leng="",
    band_str=True,
    elast_prop=True,
    optical_prop=True,
    mbj_prop=True,
    spin_orb=False,
    phonon=False,
    surf_en=False,
    def_en=False,
):
    """
    Main function to converge k-points/cut-off
    optimize structure, and run subsequent property calculations

    Args:
         mat: Poscar object with structure information
         encut: if '' then automataic convergence, else use defined fixed-cutoff
         leng: if '' then automataic convergence, else use defined fixed-line density
         band_str: if True then do band-structure calculations along high-symmetry points
         elast_prop: if True then do elastic property calculations using finite-difference
         optical_prop: if True do frequency dependent dielectric function calculations using he independent-particle (IP) approximation
         mbj_prop: if True do METAGGA-TBmBJ optical property calculations with IP approximation
         spin_orb: if True do spin-orbit calculations
         phonon: if True do phonon calculations using DFPT
         surf_en: if True do surface enrgy calculations 
         def_en: if True do defect enrgy calculations 
    Returns:
         en2: final energy
         mat_f: final structure
    """
    if encut == "":
        encut = converg_encut(encut=500, mat=mat)

    if leng == "":
        leng = converg_kpoints(length=0, mat=mat)

    kpoints = Auto_Kpoints(mat=mat, length=leng)
    isif = 2
    commen = str(mat.comment)
    lcharg = ".FALSE."
    if commen.split("@")[0] == "bulk":
        isif = 3
        lcharg = ".TRUE."
    if commen.split("@")[0] == "sbulk":
        isif = 3
    incar_dict = use_incar_dict
    incar_dict.update(
        {
            "ENCUT": encut,
            "EDIFFG": -1e-3,
            "ISIF": 3,
            "NEDOS": 5000,
            "NSW": 500,
            "NELM": 500,
            "LORBIT": 11,
            "LVTOT": ".TRUE.",
            "LVHAR": ".TRUE.",
            "ISPIN": 2,
            "LCHARG": ".TRUE.",
        }
    )
    incar = Incar.from_dict(incar_dict)
    try:
        if mat.comment.startswith("Mol"):
            incar.update({"ISIF": "2"})
    except:
        print("Mol error")
        pass
    # if commen.startswith('Surf-') :
    #   pol=check_polar(mat_f.structure)
    #   if pol==True:
    #        ase_atoms = AseAtomsAdaptor().get_atoms(mat_f.structure)
    #        COM=ase_atoms.get_center_of_mass(scaled=True)
    #        incar.update({"LDIPOL": '.TRUE.',"IDIPOL":4,"ISYM": 0,"DIPOL":COM})
    print("running smart_converge for", str(mat.comment) + str("-") + str("MAIN-RELAX"))
    cwd = str(os.getcwd())
    en2, contc = run_job(
        mat=mat,
        incar=incar,
        kpoints=kpoints,
        jobname=str("MAIN-RELAX") + str("-") + str(mat.comment),
    )
    os.chdir(cwd)
    path = str(contc.split("/CONTCAR")[0]) + str("/vasprun.xml")
    v = open(path, "r").readlines()
    for line in v:
        if "NBANDS" in line:
            nbands = int(line.split(">")[1].split("<")[0])
            print("nbands=", nbands)
            break
    strt = Structure.from_file(contc)
    mat_f = Poscar(strt)
    mat_f.comment = str(mat.comment)
    if band_str == True:
        incar_dict = use_incar_dict
        incar_dict.update(
            {
                "ISPIN": 2,
                "NEDOS": 5000,
                "LORBIT": 11,
                "IBRION": 1,
                "ENCUT": encut,
                "NBANDS": int(nbands) + 10,
            }
        )
        incar = Incar.from_dict(incar_dict)
        kpath = HighSymmKpath(mat_f.structure)
        frac_k_points, k_points_labels = kpath.get_kpoints(
            line_density=20, coords_are_cartesian=False
        )
        kpoints = Kpoints(
            comment="Non SCF run along symmetry lines",
            style=Kpoints.supported_modes.Reciprocal,
            num_kpts=len(frac_k_points),
            kpts=frac_k_points,
            labels=k_points_labels,
            kpts_weights=[1] * len(frac_k_points),
        )

        try:
            print("running MAIN-BAND")
            kpoints = mpvis.get_kpoints(mat_f.structure)
            en2B, contcB = run_job(
                mat=mat_f,
                incar=incar,
                kpoints=kpoints,
                jobname=str("MAIN-BAND") + str("-") + str(mat_f.comment),
            )
        # kpoints=mpvis.get_kpoints(mat_f.structure)
        # en2B,contcB=run_job(mat=mat_f,incar=incar,kpoints=kpoints,jobname=str('MAIN-BAND')+str('-')+str(mat_f.comment))
        except:
            print("No band str calc.")
            if str(os.getcwd) != cwd:
                print("Changing directory")
                line = str("cd ") + str(cwd)
                os.chdir(cwd)
                print(os.getcwd())

            pass
    os.chdir(cwd)
    if surf_en == True:
        incar_dict = use_incar_dict
        incar_dict.update(
            {"ENCUT": encut, "NEDOS": 5000, "IBRION": 1, "NSW": 500, "LORBIT": 11}
        )
        incar = Incar.from_dict(incar_dict)
        surf = surfer(mat=mat_f.structure, layers=3)
        for i in surf:

            try:
                print("running MAIN-BAND")
                # NSCF
                # chg_file=str(contc).replace('CONTCAR','CHGCAR')
                # print ('chrfile',chg_file)
                # shutil.copy2(chg_file,'./')
                kpoints = Auto_Kpoints(mat=i, length=leng)
                en2s, contcs = run_job(
                    mat=i,
                    incar=incar,
                    kpoints=kpoints,
                    jobname=str("Surf_en-")
                    + str(i.comment)
                    + str("-")
                    + str(mat_f.comment),
                )
                # kpoints=mpvis.get_kpoints(mat_f.structure)
                # en2B,contcB=run_job(mat=mat_f,incar=incar,kpoints=kpoints,jobname=str('MAIN-BAND')+str('-')+str(mat_f.comment))
            except:
                pass
    os.chdir(cwd)
    if def_en == True:
        incar_dict = use_incar_dict
        incar_dict.update(
            {"ENCUT": encut, "NEDOS": 5000, "IBRION": 1, "NSW": 500, "LORBIT": 11}
        )
        incar = Incar.from_dict(incar_dict)
        # surf=surfer(mat=mat_f.structure,layers=3)
        vac = vac_antisite_def_struct_gen(cellmax=3, struct=mat_f.structure)
        for i in vac:
            try:
                print("running MAIN-vac")
                kpoints = Auto_Kpoints(mat=i, length=leng)
                en2d, contcd = run_job(
                    mat=i,
                    incar=incar,
                    kpoints=kpoints,
                    jobname=str("Def_en-")
                    + str(i.comment)
                    + str("-")
                    + str(mat_f.comment),
                )
            # kpoints=mpvis.get_kpoints(mat_f.structure)
            # en2B,contcB=run_job(mat=mat_f,incar=incar,kpoints=kpoints,jobname=str('MAIN-BAND')+str('-')+str(mat_f.comment))
            except:
                pass
    os.chdir(cwd)
    # surf=surfer(mat=strt,layers=layers)
    # surf=surfer(mat=strt,layers=layers)
    if spin_orb == True:
        # chg_file=str(contc).replace('CONTCAR','CHGCAR')
        # print ('chrfile',chg_file)
        # shutil.copy2(chg_file,'./')
        incar_dict = use_incar_dict
        incar_dict.update(
            {
                "ENCUT": encut,
                "NPAR": ncores,
                "GGA_COMPAT": ".FALSE.",
                "LSORBIT": ".TRUE.",
                "IBRION": 1,
                "ISYM": 0,
                "NEDOS": 5000,
                "IBRION": 1,
                "NSW": 500,
                "LORBIT": 11,
            }
        )
        incar = Incar.from_dict(incar_dict)
        sg_mat = SpacegroupAnalyzer(mat_f.structure)
        mat_cvn = sg_mat.get_conventional_standard_structure()
        mat_cvn.sort()
        kpoints = Auto_Kpoints(mat=Poscar(mat_cvn), length=leng / 2)
        try:
            en2S, contcS = run_job(
                mat=Poscar(mat_cvn),
                incar=incar,
                kpoints=kpoints,
                jobname=str("MAIN-SOC") + str("-") + str(mat_f.comment),
            )
        except:
            pass
    os.chdir(cwd)
    if optical_prop == True:
        incar_dict = use_incar_dict
        incar_dict.update(
            {
                "NEDOS": 5000,
                "LORBIT": 11,
                "IBRION": 1,
                "ENCUT": encut,
                "NBANDS": 3 * int(nbands),
                "LOPTICS": ".TRUE.",
            }
        )
        incar = Incar.from_dict(incar_dict)
        kpoints = Auto_Kpoints(mat=mat_f, length=leng)
        try:
            en2OP, contcOP = run_job(
                mat=mat_f,
                incar=incar,
                kpoints=kpoints,
                jobname=str("MAIN-OPTICS") + str("-") + str(mat_f.comment),
            )
        except:
            pass
    os.chdir(cwd)
    if mbj_prop == True:
        incar_dict = use_incar_dict
        incar_dict.update(
            {
                "NEDOS": 5000,
                "LORBIT": 11,
                "IBRION": 1,
                "ENCUT": encut,
                "NBANDS": 3 * int(nbands),
                "LOPTICS": ".TRUE.",
                "METAGGA": "MBJ",
                "ISYM": 0,
                "SIGMA": 0.1,
            }
        )
        incar = Incar.from_dict(incar_dict)
        kpoints = Auto_Kpoints(mat=mat_f, length=leng)
        try:
            en2OP, contcOP = run_job(
                mat=mat_f,
                incar=incar,
                kpoints=kpoints,
                jobname=str("MAIN-MBJ") + str("-") + str(mat_f.comment),
            )
        except:
            pass
    os.chdir(cwd)

    if elast_prop == True:
        incar_dict = use_incar_dict
        incar_dict.update(
            {
                "NEDOS": 5000,
                "IBRION": 6,
                "ENCUT": 1.3 * float(encut),
                "ISIF": 3,
                "POTIM": 0.015,
                "NPAR": ncores,
                "ISPIN": 2,
            }
        )
        incar = Incar.from_dict(incar_dict)
        sg_mat = SpacegroupAnalyzer(mat_f.structure)
        mat_cvn = sg_mat.get_conventional_standard_structure()
        mat_cvn.sort()
        kpoints = Auto_Kpoints(mat=Poscar(mat_cvn), length=leng)
        try:
            en2E, contcE = run_job(
                mat=Poscar(mat_cvn),
                incar=incar,
                kpoints=kpoints,
                jobname=str("MAIN-ELASTIC") + str("-") + str(mat_f.comment),
            )
        except:
            pass
    os.chdir(cwd)

    if phonon == True:
        incar_dict = use_incar_dict
        incar_dict.update(
            {
                "IBRION": 8,
                "ENCUT": float(encut),
                "ISYM": 0,
                "ADDGRID": ".TRUE.",
                "EDIFF": 1e-09,
                "LORBIT": 11,
            }
        )
        incar = Incar.from_dict(incar_dict)
        kpoints = Auto_Kpoints(mat=mat_f, length=leng)
        mat_pho = make_big(poscar=mat_f, size=11.0)
        try:
            en2P, contcP = run_job(
                mat=mat_pho,
                incar=incar,
                kpoints=kpoints,
                jobname=str("MAIN-PHO8") + str("-") + str(mat_f.comment),
            )
        except:
            pass
    os.chdir(cwd)

    # if Raman_calc==True:
    #    Raman(strt=mat_f,encut=encut,length=leng)
    os.chdir(cwd)
    return en2, mat_f


def smart_vac(strt=None, tol=0.1):
    """
    Umbrell function for vacancy formation energies with convergence

    Args:
       strt: Structure object
       tol: defect energy convergence tolerance in eV
    Returns:
          def_list: list of defect energies
          def_header_list: list of defect names
           
    """
    vac_arr = []
    sg_mat = SpacegroupAnalyzer(strt)
    mat_cvn = sg_mat.get_conventional_standard_structure()
    mat_cvn.sort()

    cellmax = 1  # int(mat_cvn.composition.num_atoms)+int(mat_cvn.ntypesp)#5
    ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
    ase_atoms = ase_atoms * (cellmax, cellmax, cellmax)
    # if len(ase_atoms) >200:
    #   cellmax=1
    # else:
    #   cellmax=2
    # cellmax=int(mat_cvn.composition.num_atoms)+int(mat_cvn.ntypesp)#5
    print("type of trt is= celmmax", type(strt), cellmax)
    try:
        print("int(strt.composition.num_atoms)", int(strt.composition.num_atoms))
        print(int(strt.ntypesp))
    except:
        pass
    # cellmax=int(strt.composition.num_atoms)+int(strt.ntypesp)
    vac_done = 0
    vac = vac_antisite_def_struct_gen(cellmax=cellmax, struct=strt)
    def_list = [100000 for y in range(len(vac) - 1)]
    while vac_done != 1:
        vac = vac_antisite_def_struct_gen(cellmax=cellmax, struct=strt)
        if vac not in vac_arr:
            vac_arr.append(vac)
            print(
                "in smart_vac(strt=None), cellmax,vac,vac_done=", cellmax, vac, vac_done
            )
            def_list2, header_list = def_energy(vac=vac)
            diff = matrix(def_list) - matrix(def_list2)
            diff_arr = np.array(diff).flatten()
            print("in smart_vac(strt=None diff_arr=", diff_arr)
            if any(diff_arr) > tol:
                # for el in diff_arr:
                #     if abs(el)>tol :
                # print ("in smart_vac(strt=None abs_el=",abs(el))
                vac_done = 0
                cellmax = cellmax + 1
                ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
                ase_atoms = ase_atoms * (cellmax, cellmax, cellmax)
                if len(ase_atoms) > 100:
                    vac_done = 1
                def_list = def_list2
            else:
                vac_done = 1
    #        cellmax=cellmax+1
    return def_list, header_list


def def_energy(vac=[]):
    """
    Calculation of vacancy formation energies

    Args:
        vac: list of Poscar vacancy structure objects
    Returns:
          def_list: list of defect energies
          def_header_list: list of defect names
    """

    def_list = []
    header_list = []
    fi = str("DEF.INFO")
    f = open(fi, "w")
    for v in vac:
        enp, contc = smart_converge(mat=v)
        # enp,contc=run_job(mat=v,incar=incar,kpoints=kpoints,jobname=str(v.comment))
        strt = Structure.from_file(contc)
        comm = str(v.comment)
        print("running def_energy for =", comm)
        header_list.append(comm)
        if comm.split("@")[0] == "bulk":
            gs_energy = float(enp) * float(strt.composition.num_atoms)
            print("in def_energy gs_energy for", comm, gs_energy)
        else:
            chem_pot = sum_chem_pot(strt)
            if comm.startswith("intl"):
                chem_pot = 0.0 - float(chem_pot)
            def_en = (
                (float(enp) * float(strt.composition.num_atoms))
                - float(gs_energy)
                + chem_pot
            )
            print("in def_energy def_en for", comm, def_en, "chem_pot", chem_pot)
            def_list.append(def_en)
            print(v.comment, "=", def_en)
            line = str(v.comment) + str("=") + str(def_en) + "\n"
            f.write(line)
    f.close()
    return def_list, header_list


def smart_surf(strt=None, tol=0.1):
    """
    Umbrell function for surface energies with convergence

    Args:
       strt: Structure object
       tol: surface energy convergence tolerance in eV
    Returns:
          surf_list: list of surface energies
          surf_header_list: list of surface names
           
    """
    sg_mat = SpacegroupAnalyzer(strt)
    mat_cvn = sg_mat.get_conventional_standard_structure()
    mat_cvn.sort()
    layers = 2
    indices = get_symmetrically_distinct_miller_indices(mat_cvn, 1)
    ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
    for i in indices:
        ase_slab = surface(ase_atoms, i, layers)
        ase_slab.center(vacuum=15, axis=2)
        if len(ase_slab) < 50:
            layers = 3
    surf_arr = []
    surf_done = 0
    surf = surfer(mat=strt, layers=layers)
    surf_list = [100000 for y in range(len(surf) - 1)]
    print("in smart_surf :surf,surf_list=", surf, surf_list)
    while surf_done != 1:
        layers = layers + 1
        indices = get_symmetrically_distinct_miller_indices(mat_cvn, 1)
        ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
        for i in indices:
            ase_slab = surface(ase_atoms, i, layers)
            ase_slab.center(vacuum=15, axis=2)
            if len(ase_slab) > 100:
                surf_done = 1
            if (ase_slab.get_cell()[2][2]) > 40:
                surf_done = 1
        surf = surfer(mat=strt, layers=layers)
        if surf not in surf_arr:
            surf_arr.append(surf)
            surf_list2, surf_header_list = surf_energy(surf=surf)
            print("in smart_surf :surf2,surf_list2=", surf_list2, surf_header_list)
            diff = matrix(surf_list) - matrix(surf_list2)
            print(
                "in smart_surf :surf3,surf_list3=",
                matrix(surf_list),
                matrix(surf_list2),
            )
            diff_arr = np.array(diff).flatten()
            if any(diff_arr) > tol:
                # for el in diff_arr:
                #    if abs(el)>tol :
                #        print ("in smart_surf :abs el=",abs(el))
                surf_done = 0
                surf_list = surf_list2
            else:
                surf_done = 1
    return surf_list, surf_header_list


def surf_energy(surf=[]):
    """
    Helper function for surface energies

    Args:
        surf: list of Poscar surface objects
    Returns:
          surf_list: list of surface energies
          surf_header_list: list of surface names
          
    """

    surf_list = []
    surf_header_list = []
    fi = str("SURF.INFO")
    f = open(fi, "w")
    for s in surf:
        enp, contc = smart_converge(mat=s)
        strt = Structure.from_file(contc)
        m = strt.lattice.matrix
        if s.comment.split("@")[0] == "sbulk":

            gs_energy = enp
            print("in surf_energy gs_energy for", s.comment, gs_energy)
        else:
            surf_en = (
                16
                * (
                    float(enp) * (strt.composition.num_atoms)
                    - float(strt.composition.num_atoms) * float(gs_energy)
                )
                / (2 * np.linalg.norm(np.cross(m[0], m[1])))
            )
            print("in surf_energy surf_en for", s.comment, surf_en)
            surf_list.append(surf_en)
            print(s.comment, "=", surf_en)
            line = str(s.comment) + str("=") + str(surf_en) + "\n"
            f.write(line)
    f.close()
    return surf_list, surf_header_list


def make_big(poscar=None, size=11.0):
    """
    Helper function to make supercell
 
    Args:
      poscar: Poscar object
      size: simulation size in Angstrom
    Returns:
       big: Poscar supercell object
    """
    struct = poscar.structure
    comm = poscar.comment
    a, b, c = struct.lattice.abc
    struct.make_supercell(
        [
            int(float(size) / float(a)) + 1,
            int(float(size) / float(b)) + 1,
            int(float(size) / float(c)) + 1,
        ]
    )
    big = Poscar(struct)
    big.comment = str(comm)
    return big


def main_func(mpid="", jid="", mat=None, enforc_cvn=False):
    """
    Main function to carry out property calculations

    Args:
        mpid: materialsproject id
        jid: jarvis-dft id
        mat: Poscar object
        enforc_cvn:  whether or not enforce conventional cell input structure
    Returns:
       en: final energy
       final: final structure
    """

    if mpid != "" or jid != "":
        data = loadfn(json_dat, cls=MontyDecoder)
        for d in data:
            mmpid = str(d["mpid"])
            jjid = str(d["mpid"])
            if mmpid == mpid:
                fin = d["structure"]
                break
            if jjpid == jid:
                fin = d["structure"]
                break

        strt = fin
        sg_mat = SpacegroupAnalyzer(strt)
        mat_cvn = sg_mat.get_conventional_standard_structure()
        mat_cvn.sort()
        if (
            int(strt.composition._natoms) == int(mat_cvn.composition._natoms)
            and enforc_cvn == True
        ):
            mat = Poscar(mat_cvn)
        else:
            mat = Poscar(strt)

        mpid = mpid.replace("-", "_")
        mpid = str("bulk@") + str(mpid)
        mat.comment = mpid

    en, final = smart_converge(mat=mat)
    print(en, final)


# if __name__ == '__main__':
# struct=Structure.from_file("POSCAR")
# pos=Poscar.from_file('POSCAR')
# #pos.comment='Surf-mp-C'
# main_func(mat=pos)
