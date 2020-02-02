from jarvis.io.vasp.outputs import Outcar
from jarvis.io.vasp.inputs import Poscar, Incar, Potcar
from jarvis.db.jsonutils import loadjson, dumpjson
import os
import shutil
from jarvis.core.kpoints import Kpoints3D as Kpoints


class VaspJobs(object):
    def __init__(
        self,
        poscar=None,
        kpoints=None,
        incar=None,
        potcar = None,
        pot_type='POT_GGA_PAW_PBE',
        ncores=16,
        ncpu=1,
        jobname="",
        copy_files=[],
    ):
        self.poscar = poscar
        self.kpoints = kpoints
        self.incar = incar
        self.potcar=potcar
        self.pot_type=pot_type
        self.ncores = ncores
        self.ncpu = ncpu
        self.copy_files = copy_files
        self.jobname = jobname

    def main(self):
        #poscar=self.poscar
        #incar=self.incar
        #kpoints=self.kpoints
        #copy_files=self.copy_files


        cwd = str(os.getcwd())
        if self.jobname =="":
              jobname = str(self.poscar.comment)
        job_dir = str(self.jobname)
        run_file = str(os.getcwd()) + str("/") + str(self.jobname) + str(".json")
        run_dir = str(os.getcwd()) + str("/") + str(self.jobname)
        if self.poscar.comment.startswith("Surf"):
            [a, b, c] = self.kpoints.kpts[0]
            self.kpoints.kpts = [[a, b, 1]]
            try:
                pol = self.poscar.atoms.check_polar
                if pol == True:
                    COM = self.poscar.atoms.get_center_of_mass()
                    print("COM=", COM)
                    print("Found polar surface, will be setting dipole corrections")
                    self.incar.update(
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
                    print("Polar surface encountered in run_job", poscar.comment)
            except:
                pass
        wait = False
        json_file = str(self.jobname) + str(".json")
        print("json should be here=", str(os.getcwd()) + str("/") + str(json_file))
        # print ('json should be=',json_file,run_file,os.getcwd())
        if os.path.exists(str(os.getcwd()) + str("/") + str(json_file)):
            try:
                data_cal = loadjson(str(os.getcwd()) + str("/") + str(json_file))
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
            #print("Setting up POTCAR")
            if self.potcar is None:
              new_symb = list(set(self.poscar.atoms.elements))
              self.potcar = Potcar(elements=new_symb, pot_type=self.pot_type)
            if not os.path.exists(run_dir):
                print("Starting new job")
                os.makedirs(run_dir)
                os.chdir(run_dir)
                self.incar.write_file("INCAR")
                self.potcar.write_file("POTCAR")
                self.kpoints.write_file("KPOINTS")
                self.poscar.write_file("POSCAR")
                for i in self.copy_files:
                    print("copying", i)
                    shutil.copy2(i, "./")
                # f=open('job.out','w')

                cmd = (
                    str("/users/knc6/VASP/vasp54/src/vasp.5.4.1DobbySOC2/bin/vasp_ncl")
                    + "\n"
                )

                os.system(cmd)
                if os.path.isfile("OUTCAR"):
                    try:
                        wait = Outcar(
                            "OUTCAR"
                        ).converged  # Vasprun("vasprun.xml").converged
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
                    self.incar.write_file("INCAR")
                    self.kpoints.write_file("KPOINTS")
                    self.poscar.write_file("POSCAR")
                    try:
                        if self.potcar is None:
                             new_symb = list(set(self.poscar.atoms.elements))
                             self.potcar = Potcar(elements=new_symb, pot_type=self.pot_type)
                        self.potcar.write_file("POTCAR")
                        print("FOUND OLD CONTCAR in", os.getcwd())
                        # old_contcar = Poscar.from_file("CONTCAR")
                        # old_contcar.write_file('POSCAR')
                        copy_cmd = str("cp CONTCAR POSCAR")
                        print("copy_cmd=", copy_cmd)
                        if "ELAST" not in jobname and "LEPSILON" not in jobname:
                            # Because in ELASTIC calculations structures are deformed
                            os.system(copy_cmd)
                        # time.sleep(3)
                    except:
                        pass
                    for i in self.copy_files:
                        print("copying", i)
                        shutil.copy2(i, "./")

                    cmd = (
                        str(
                            "/users/knc6/VASP/vasp54/src/vasp.5.4.1DobbySOC2/bin/vasp_ncl"
                        )
                        + "\n"
                    )
                    os.system(cmd)
                    if os.path.isfile("OUTCAR"):
                        try:
                            wait = Outcar(
                                "OUTCAR"
                            ).converged  # Vasprun("vasprun.xml").converged
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
                    "poscar_initial": poscar.atoms.as_dict(),
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


if __name__ == "__main__":

    p = Poscar.from_file(
        "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-RELAX-bulk@mp_541837/CONTCAR"
    )
    #p.write_file('POSCAR')
    inc = Incar.from_file(
        "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-RELAX-bulk@mp_541837/INCAR"
    )
    #inc.write_file('INCAR')
    kp = Kpoints().automatic_length_mesh(lattice_mat=p.atoms.lattice_mat)
    #kp.write_file('KPOINTS')
    # pot=Potcar(elements=p.atoms.elements)
    new_symb = list(set(p.atoms.elements))
    potcar = Potcar(elements=new_symb)
    job = VaspJobs(poscar=p, kpoints=kp, incar=inc, jobname="testt").main()
