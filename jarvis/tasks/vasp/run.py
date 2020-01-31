from jarvis.io.vasp.outputs import Outcar
from jarvis.io.vasp.inputs import Poscar,Incar,Potcar
from jarvis.db.jsonutils import loadjson,dumpjson
import os
import shutil
from jarvis.core.kpoints import Kpoints3D as Kpoints

class VaspJobs(object):
   def __init__(self,atoms=None, kpoints=None, incar=None, pot_type='POT_GGA_PAW_PBE', ncores=16, ncpu=1, jobname='',copy_files=[]):
       self.atoms = atoms
       self.kpoints=kpoints
       self.incar=incar
       self.potcar=potcar
       self.ncores=ncores
       self.ncpu=ncpu

   def main(self):
     

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
            pol = mat.atoms.check_polar 
            if pol == True:
                COM = mat.get_center_of_mass()
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
            data_cal = loadjson(
                str(os.getcwd()) + str("/") + str(json_file)
            )
            tmp_outcar = (
                str(os.getcwd())
                + str("/")
                + str(json_file.split(".json")[0])
                + str("/OUTCAR")
            )
            print("outcar is", tmp_outcar)
            wait = Outcar(tmp_outcar).converged  # True
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
        new_symb = mat.elements
        try:
            potcar = Potcar(symbols=new_symb, pot_type=pot_type)
        except:
            print(
                "JARVIS-ERROR: Could not set POTCAR, check JARVIS_VASP_PSP_DIR"
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


            os.system(cmd)
            if os.path.isfile("OUTCAR"):
                try:
                    wait = Outcar("OUTCAR").converged  # Vasprun("vasprun.xml").converged
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
                    #old_contcar = Poscar.from_file("CONTCAR")
                    # old_contcar.write_file('POSCAR')
                    copy_cmd = str("cp CONTCAR POSCAR")
                    print("copy_cmd=", copy_cmd)
                    if "ELAST" not in jobname and 'LEPSILON' not in jobname:
                        # Because in ELASTIC calculations structures are deformed
                        os.system(copy_cmd)
                    # time.sleep(3)
                except:
                    pass
                for i in copy_file:
                    print("copying", i)
                    shutil.copy2(i, "./")

                cmd = str("python  first_cust.py >out_dat") + "\n"
                os.system(cmd)
                if os.path.isfile("OUTCAR"):
                    try:
                        wait = Outcar("OUTCAR").converged  # Vasprun("vasprun.xml").converged
                        # wait=Vasprun("vasprun.xml").converged
                    except:
                        pass
    f_energy = "na"
    enp = "na"
    contcar = str(os.getcwd()) + str("/") + str("CONTCAR")
    final_str = Poscar.from_file(contcar).atoms
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
