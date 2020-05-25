"""
Class for running VASP jobs
"""

from jarvis.io.vasp.outputs import Outcar, Vasprun
from jarvis.io.vasp.inputs import Poscar, Incar, Potcar
from jarvis.db.jsonutils import loadjson, dumpjson
from jarvis.core.kpoints import Kpoints3D as Kpoints
from jarvis.tasks.queue_jobs import Queue
from jarvis.analysis.structure.spacegroup import Spacegroup3D
import subprocess
import json
import os
import shutil


class JobFactory(object):
    def __init__(self, name="", use_incar_dict={}, pot_type=None):
        """
        Generic class for running variations of VASP calculations
        
        Args:
        
            name : generic name
            
            use_incar_dict : dictionary with INCAR parameters that would be repreated
            
            pot_type : pseudopotential type
        """
        self.name = name
        self.use_incar_dict = use_incar_dict
        self.pot_type = pot_type

    def all_optb88vdw_props(self, mat=None):
        """
        Used for OptB88vdW functional based high-throughput calculations
        This will converge k-points, cut-offs, and then carry several property calculations.
        
        Args:
        
            mat : Poscar object
        """
        optb88 = GenericIncars().optb88vdw()
        job = JobFactory(use_incar_dict=optb88.incar, pot_type=optb88.pot_type)
        encut = job.converg_encut(mat=mat)
        length = job.converg_kpoint(mat=mat)
        energy, contcar_path = job.optimize_geometry(
            mat=mat, encut=encut, length=length
        )
        optimized_mat = Poscar.from_file(contcar_path)
        vrun = Vasprun(contcar_path.replace("CONTCAR", "vasprun.xml"))
        chg_path = contcar_path.replace("CONTCAR", "CHGCAR")
        nbands = int(vrun.all_input_parameters["NBANDS"])
        enB, contcB = job.band_structure(
            mat=optimized_mat,
            encut=encut,
            line_density=20,
            nbands=2 * nbands,
            copy_prev_chgcar=chg_path,
        )
        enL, contcL = job.loptics(
            mat=optimized_mat, encut=encut, nbands=2 * nbands, length=length
        )
        enM, contcM = job.mbj_loptics(
            mat=optimized_mat, encut=encut, nbands=2 * nbands, length=length
        )
        enE, contcE = job.elastic(
            mat=optimized_mat, encut=encut, nbands=2 * nbands, length=length
        )

    def elastic(
        self, mat=None, encut=None, nbands=None, potim=0.015, npar=None, length=20
    ):
        """
        Used for elastic property calculations using IBRION = 6
        Enforces conventional standard structure
        
        Args:
        
            mat :  Poscar object
            
            encut :  Plane-wave cut-off, 1.3 times will be used
            
            nbands : number of bands, generally high-value recommended
            
            npar : NPAR tag, see VASP manual, set it as number of cores
            
            length :  K-points in length unit
        """

        incar = self.use_incar_dict
        cvn = Spacegroup3D(mat.atoms).conventional_standard_structure
        comment = mat.comment
        p = Poscar(cvn, comment=comment)

        if npar is not None:
            incar.update({"NPAR": npar})

        if nbands is not None:
            nbands = int(nbands * 3)
            incar.update({"NBANDS": nbands})
        data = {
            "ENCUT": 1.3 * float(encut),
            "NEDOS": 5000,
            "ISIF": 3,
            "POTIM": potim,
            "ISPIN": 2,
            "IBRION": 6,
            "LCHARG": ".FALSE.",
        }
        incar.update(data)
        kpoints = Kpoints().automatic_length_mesh(
            lattice_mat=p.atoms.lattice_mat, length=length
        )  # Auto_Kpoints(mat=mat, length=length)

        en, contcar = VaspJobs(
            poscar=p,
            incar=incar,
            pot_type=self.pot_type,
            kpoints=kpoints,
            jobname=str("MAIN-ELASTIC-") + str(p.comment.split()[0]),
        ).runjob()

        return en, contcar

    def mbj_loptics(self, mat=None, encut=None, nbands=None, length=20):
        """
        Used for TBmBJ meta-GGA calculation
        
        Args:
        
            mat :  Poscar object
            
            encut :  Plane-wave cut-off, 1.3 times will be used
            
            nbands : number of bands, increased to threee times
            
            length :  K-points in length unit
        """
        incar = self.use_incar_dict

        if nbands is not None:
            nbands = int(nbands * 3)
            incar.update({"NBANDS": nbands})
        data = {
            "ENCUT": encut,
            "NEDOS": 5000,
            "NELM": 500,
            "LORBIT": 11,
            "ISPIN": 2,
            "METAGGA": "MBJ",
            "SIGMA": 0.1,
            "ISYM": 0,
            "LOPTICS": ".TRUE.",
            "IBRION": 1,
            "LCHARG": ".FALSE.",
        }
        incar.update(data)
        kpoints = Kpoints().automatic_length_mesh(
            lattice_mat=mat.atoms.lattice_mat, length=length
        )  # Auto_Kpoints(mat=mat, length=length)

        en, contcar = VaspJobs(
            poscar=mat,
            incar=incar,
            pot_type=self.pot_type,
            kpoints=kpoints,
            jobname=str("MAIN-MBJ-") + str(mat.comment.split()[0]),
        ).runjob()

        return en, contcar

    def loptics(self, mat=None, encut=None, nbands=None, length=20):
        """
        Used in linear-optics calculations
        
        Args:
        
            mat :  Poscar object
            
            encut :  Plane-wave cut-off, 1.3 times will be used
            
            nbands : number of bands, increased to threee times
            
            length :  K-points in length unit
        """
        incar = self.use_incar_dict

        if nbands is not None:
            nbands = int(nbands * 3)
            incar.update({"NBANDS": nbands})
        data = {
            "ENCUT": encut,
            "NEDOS": 5000,
            "NELM": 500,
            "LORBIT": 11,
            "ISPIN": 2,
            "LOPTICS": ".TRUE.",
            "IBRION": 1,
            "LCHARG": ".FALSE.",
        }
        incar.update(data)
        kpoints = Kpoints().automatic_length_mesh(
            lattice_mat=mat.atoms.lattice_mat, length=length
        )  # Auto_Kpoints(mat=mat, length=length)

        en, contcar = VaspJobs(
            poscar=mat,
            incar=incar,
            pot_type=self.pot_type,
            kpoints=kpoints,
            jobname=str("MAIN-OPTICS-") + str(mat.comment.split()[0]),
        ).runjob()

        return en, contcar

    def band_structure(
        self, mat=None, encut=None, line_density=20, nbands=None, copy_prev_chgcar=None
    ):
        """
        Used in band-structure calculations
        
        Args:
        
            mat :  Poscar object
            
            encut :  Plane-wave cut-off, 1.3 times will be used
            
            nbands : number of bands, increased to threee times
            
            line_density :  number of k-points between two high-symmetry k-points
            
            copy_prev_chgcar :  path of CHGCAR file for Non-SCF step
        """
        incar = self.use_incar_dict
        if copy_prev_chgcar is not None:
            shutil.copy2(copy_prev_chgcar, ".")

        if nbands is not None:
            nbands = int(nbands * 1.3)
            incar.update({"NBANDS": nbands})
        data = {
            "ENCUT": encut,
            "NEDOS": 5000,
            "NELM": 500,
            "LORBIT": 11,
            "ISPIN": 2,
            "IBRION": 1,
            "LCHARG": ".FALSE.",
        }
        incar.update(data)
        kpoints = Kpoints().kpath(mat.atoms, line_density=line_density)

        en, contcar = VaspJobs(
            poscar=mat,
            incar=incar,
            pot_type=self.pot_type,
            kpoints=kpoints,
            jobname=str("MAIN-BAND-") + str(mat.comment.split()[0]),
        ).runjob()

        return en, contcar

    def optimize_geometry(self, mat=None, encut=None, length=None):
        """
        Used in optimizing lattice-parameter and internal psotions
        
        Args:
        
            mat :  Poscar object
            
            encut :  Plane-wave cut-off
            
            length :  K-points in length unit
        """
        incar = self.use_incar_dict
        data = {
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
        incar.update(data)
        kpoints = Kpoints().automatic_length_mesh(
            lattice_mat=mat.atoms.lattice_mat, length=length
        )  # Auto_Kpoints(mat=mat, length=length)
        en, contcar = VaspJobs(
            poscar=mat,
            incar=incar,
            pot_type=self.pot_type,
            kpoints=kpoints,
            jobname=str("MAIN-RELAX-") + str(mat.comment),
        ).runjob()
        return en, contcar

    def converg_encut(self, encut=500, mat=None, starting_length=10, tol=0.001):
        """
        Function to converg plane-wave cut-off
        
        Args:
        
            encut: intial cutoff
            
            mat: Poscar object
            
        Returns:
        
               encut: converged cut-off
        """
        pot_type = self.pot_type
        en1 = -10000
        encut1 = encut
        convg_encut1 = False
        convg_encut2 = False

        while convg_encut2 != True:
            # while convg_encut1 !=True and  convg_encut2 !=True:
            # tol = 0.001  # change 0.001
            encut_list = []
            encut_list.append(encut)
            length = starting_length
            encut1 = encut + 50
            incar_dict = self.use_incar_dict
            # print ('incar_dict',incar_dict)
            incar_dict.update({"ENCUT": encut})
            # print (use_incar_dict)
            incar = incar_dict
            kpoints = Kpoints().automatic_length_mesh(
                lattice_mat=mat.atoms.lattice_mat, length=length
            )  # Auto_Kpoints(mat=mat, length=length)
            print(
                "running smart_converge for",
                str(mat.comment) + str("-") + str("ENCUT") + str("-") + str(encut),
            )
            en2, contc = VaspJobs(
                poscar=mat,
                incar=incar,
                pot_type=pot_type,
                kpoints=kpoints,
                jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut),
            ).runjob()
            while abs(en2 - en1) > tol:
                en1 = en2
                encut1 = encut + 50
                encut_list.append(encut)
                print("Incrementing encut", encut)
                # incar_dict = self.use_incar_dict
                incar_dict.update({"ENCUT": encut1})
                # incar_dict["ENCUT"]= encut1
                # incar = Incar.from_dict(incar_dict)
                print(
                    "running smart_converge for",
                    str(mat.comment) + str("-") + str("ENCUT") + str("-") + str(encut),
                )
                en2, contc = VaspJobs(
                    poscar=mat,
                    incar=incar,
                    pot_type=pot_type,
                    kpoints=kpoints,
                    jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut),
                ).runjob()
            convg_encut1 = True

            # Some extra points to check
            print("Some extra points to check for ENCUT")

            encut2 = encut1 + 50
            incar.update({"ENCUT": encut2})
            # incar_dict["ENCUT"]= encut2
            en3, contc = VaspJobs(
                poscar=mat,
                incar=incar,
                pot_type=pot_type,
                kpoints=kpoints,
                jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut2),
            ).runjob()

            encut3 = encut2 + 50
            # incar["ENCUT"] = encut3
            incar.update({"ENCUT": encut3})
            # incar_dict["ENCUT"]= encut3
            en4, contc = VaspJobs(
                poscar=mat,
                incar=incar,
                pot_type=pot_type,
                kpoints=kpoints,
                jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut3),
            ).runjob()

            encut4 = encut3 + 50
            incar.update({"ENCUT": encut4})
            # incar_dict["ENCUT"]= encut4
            en5, contc = VaspJobs(
                poscar=mat,
                incar=incar,
                pot_type=pot_type,
                kpoints=kpoints,
                jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut4),
            ).runjob()

            encut5 = encut4 + 50
            # incar["ENCUT"] = encut5
            incar.update({"ENCUT": encut5})
            # incar_dict["ENCUT"]= encut5
            en6, contc = VaspJobs(
                poscar=mat,
                pot_type=pot_type,
                incar=incar,
                kpoints=kpoints,
                jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut5),
            ).runjob()

            encut6 = encut5 + 50
            # incar["ENCUT"] = encut6
            incar.update({"ENCUT": encut6})
            # incar_dict["ENCUT"]= encut6
            en7, contc = VaspJobs(
                poscar=mat,
                pot_type=pot_type,
                incar=incar,
                kpoints=kpoints,
                jobname=str("ENCUT") + str(mat.comment) + str("-") + str(encut6),
            ).runjob()

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

    def converg_kpoint(self, length=0, mat=None, encut=500, tol=0.001):
        """
        Function to converg K-points
        
        Args:
        
            lenght: K-point line density
            
            mat: Poscar object with structure information
            
        Returns:
        
               length1: K-point line density
        """

        pot_type = self.pot_type
        en1 = -10000
        convg_kp1 = False
        convg_kp2 = False
        length1 = length
        kp_list = []
        while convg_kp2 != True:
            # while convg_kp1 !=True and  convg_kp2 !=True:
            incar = self.use_incar_dict
            incar.update({"ENCUT": encut})
            # incar_dict["ENCUT"]= encut
            length1 = length1 + 5
            print("Incrementing length", length1)
            kpoints = Kpoints().automatic_length_mesh(
                lattice_mat=mat.atoms.lattice_mat, length=length
            )  # Auto_Kpoints(mat=mat, length=length)
            mesh = kpoints.kpts[0]
            if mesh not in kp_list:
                kp_list.append(mesh)
                en2, contc = VaspJobs(
                    poscar=mat,
                    incar=incar,
                    pot_type=pot_type,
                    kpoints=kpoints,
                    jobname=str("KPOINTS") + str(mat.comment) + str("-") + str(length1),
                ).runjob()

                while abs(en2 - en1) > tol:
                    en1 = en2
                    print("Incrementing length", length1)
                    while mesh in kp_list:
                        length1 = length1 + 5
                        ##Assuming you are not super unlucky
                        # kpoints = Auto_Kpoints(mat=mat, length=length1)
                        kpoints = Kpoints().automatic_length_mesh(
                            lattice_mat=mat.atoms.lattice_mat, length=length1
                        )  # Auto_Kpoints(mat=mat, length=length)
                        mesh = kpoints.kpts[0]

                    kpoints = Kpoints().automatic_length_mesh(
                        lattice_mat=mat.atoms.lattice_mat, length=length1
                    )  # Auto_Kpoints(mat=mat, length=length)
                    mesh = kpoints.kpts[0]
                    if mesh not in kp_list:
                        kp_list.append(mesh)
                        en2, contc = VaspJobs(
                            poscar=mat,
                            incar=incar,
                            pot_type=pot_type,
                            kpoints=kpoints,
                            jobname=str("KPOINTS")
                            + str(mat.comment)
                            + str("-")
                            + str(length1),
                        ).runjob()
                    else:
                        length1 = length1 + 5
                        ##Assuming you are not super unlucky
                        #            kpoints = Auto_Kpoints(mat=mat, length=length1)
                        kpoints = Kpoints().automatic_length_mesh(
                            lattice_mat=mat.atoms.lattice_mat, length=length1
                        )  # Auto_Kpoints(mat=mat, length=length)
                        mesh = kpoints.kpts[0]
                        kp_list.append(mesh)
                        en2, contc = Vaspjobs(
                            mat=mat,
                            incar=incar,
                            kpoints=kpoints,
                            jobname=str("KPOINTS")
                            + str(mat.comment)
                            + str("-")
                            + str(length1),
                        ).runjob()
                convg_kp1 = True

                # Some extra points to check
                print("Some extra points to check for KPOINTS")
                length3 = length1 + 5
                # kpoints = Auto_Kpoints(mat=mat, length=length3)
                kpoints = Kpoints().automatic_length_mesh(
                    lattice_mat=mat.atoms.lattice_mat, length=length3
                )  # Auto_Kpoints(mat=mat, length=length)
                mesh = kpoints.kpts[0]
                kp_list.append(mesh)
                en3, contc = VaspJobs(
                    poscar=mat,
                    pot_type=pot_type,
                    incar=incar,
                    kpoints=kpoints,
                    jobname=str("KPOINTS") + str(mat.comment) + str("-") + str(length3),
                ).runjob()

                length4 = length3 + 5
                # kpoints = Auto_Kpoints(mat=mat, length=length4)
                kpoints = Kpoints().automatic_length_mesh(
                    lattice_mat=mat.atoms.lattice_mat, length=length4
                )  # Auto_Kpoints(mat=mat, length=length)
                mesh = kpoints.kpts[0]
                kp_list.append(mesh)
                en4, contc = VaspJobs(
                    poscar=mat,
                    pot_type=pot_type,
                    incar=incar,
                    kpoints=kpoints,
                    jobname=str("KPOINTS") + str(mat.comment) + str("-") + str(length4),
                ).runjob()

                length5 = length4 + 5
                # kpoints = Auto_Kpoints(mat=mat, length=length5)
                kpoints = Kpoints().automatic_length_mesh(
                    lattice_mat=mat.atoms.lattice_mat, length=length5
                )  # Auto_Kpoints(mat=mat, length=length)
                mesh = kpoints.kpts[0]
                kp_list.append(mesh)
                en5, contc = VaspJobs(
                    poscar=mat,
                    incar=incar,
                    pot_type=pot_type,
                    kpoints=kpoints,
                    jobname=str("KPOINTS") + str(mat.comment) + str("-") + str(length5),
                ).runjob()

                length6 = length5 + 5
                # kpoints = Auto_Kpoints(mat=mat, length=length6)
                kpoints = Kpoints().automatic_length_mesh(
                    lattice_mat=mat.atoms.lattice_mat, length=length6
                )  # Auto_Kpoints(mat=mat, length=length)
                mesh = kpoints.kpts[0]
                kp_list.append(mesh)
                en6, contc = VaspJobs(
                    poscar=mat,
                    incar=incar,
                    pot_type=pot_type,
                    kpoints=kpoints,
                    jobname=str("KPOINTS") + str(mat.comment) + str("-") + str(length6),
                ).runjob()
                length7 = length6 + 5
                # kpoints = Auto_Kpoints(mat=mat, length=length7)
                kpoints = Kpoints().automatic_length_mesh(
                    lattice_mat=mat.atoms.lattice_mat, length=length7
                )  # Auto_Kpoints(mat=mat, length=length)
                mesh = kpoints.kpts[0]
                kp_list.append(mesh)
                en7, contc = VaspJobs(
                    poscar=mat,
                    incar=incar,
                    pot_type=pot_type,
                    kpoints=kpoints,
                    jobname=str("KPOINTS") + str(mat.comment) + str("-") + str(length7),
                ).runjob()

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


class VaspJobs(object):
    def __init__(
        self,
        poscar=None,
        kpoints=None,
        incar=None,
        potcar=None,
        vasp_cmd="mpirun /users/knc6/VASP/vasp54/src/vasp.5.4.1DobbySOC2/bin/vasp_std",
        output_file="vasp.out",
        stderr_file="std_err.txt",
        jobname="test",
        pot_type=None,
        copy_files=["/users/knc6/bin/vdw_kernel.bindat"],
        attempts=5,
    ):
        """
        Class defninig a typical VASP calculation
        
        Args:
        
            poscar :  Poscar object
            
            incar : Incar object
            
            kpoints : Kpoints object
            
            potcar : Potcar object
            
            vasp_cmd :  path to vasp executable
            
            output_file : standard output file
            
            stderr_file : standard error output file
            
            jobname : job name
            
            pot_type :  pseudopotential type
            
            copy_files :  file(s) to be copied
            
            attempts :  used in error handling
            
        """
        self.poscar = poscar
        self.kpoints = kpoints
        self.incar = incar
        self.potcar = potcar
        self.pot_type = pot_type
        self.vasp_cmd = vasp_cmd
        self.copy_files = copy_files
        self.attempts = attempts
        self.output_file = output_file
        self.stderr_file = stderr_file
        self.jobname = jobname
        if self.potcar is None:
            if self.pot_type is None:
                ValueError("Either pass the Potcar object or provide pot_type")

            new_symb = []
            for i in self.poscar.atoms.elements:
                if i not in new_symb:
                    new_symb.append(i)
            self.potcar = Potcar(elements=new_symb, pot_type=self.pot_type)

    def run(self):
        """
        Use subprocess to tun a job
        """
        with open(self.output_file, "w") as f_std, open(
            self.stderr_file, "w", buffering=1
        ) as f_err:
            # use line buffering for stderr
            p = subprocess.Popen(self.vasp_cmd, shell=True, stdout=f_std, stderr=f_err)
            p.wait()
        return p

    def write_jobsub_py(self, filename="jobsub.py"):
        """
        Writes a generic python file for running jobs
        """
        f = open(filename, "w")
        f.write("%s\n" % "from jarvis.io.vasp.inputs import Poscar, Incar, Potcar")
        f.write("%s\n" % "from jarvis.core.kpoints import Kpoints3D as Kpoints")
        f.write("%s\n" % 'pos=Poscar.from_file("POSCAR")')
        f.write("%s\n" % 'inc=Poscar.from_file("INCAR")')
        f.write("%s\n" % 'pot=Potcar.from_file("POTCAR")')
        f.write("%s\n" % 'kp=Kpoints.from_file("KPOINTS")')
        line = (
            str("job = VaspJobs(poscar=pos, kpoints=kp,potcar=pot, incar=inc, jobname=")
            + str(self.jobname)
            + str(").runjob()")
        )
        f.write("%s\n" % line)
        f.close()

    def runjob(self):
        """
        Main function for running a generic VASP calculation
        """

        # poscar=self.poscar
        # incar=self.incar
        # kpoints=self.kpoints
        # copy_files=self.copy_files

        cwd = str(os.getcwd())
        if self.jobname == "":
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
        attempt = 0
        while wait == False:
            attempt = attempt + 1
            if attempt == self.attempts:
                wait = True
            # print("Setting up POTCAR")
            # if self.potcar is None:
            #  new_symb = list(set(self.poscar.atoms.elements))
            #  self.potcar = Potcar(elements=new_symb, pot_type=self.pot_type)
            if not os.path.exists(run_dir):
                print("Starting new job")
                os.makedirs(run_dir)
                os.chdir(run_dir)
                self.poscar.write_file("POSCAR")
            else:
                os.chdir(run_dir)
                if os.path.isfile("OUTCAR"):
                    try:
                        wait = main_outcar("OUTCAR")  # Vasprun("vasprun.xml").converged
                        # wait=Vasprun("vasprun.xml").converged
                    except:
                        pass
                    try:
                        self.potcar.write_file("POTCAR")
                        print("FOUND OLD CONTCAR in", os.getcwd())
                        copy_cmd = str("cp CONTCAR POSCAR")
                        self.poscar.write_file("POSCAR")
                        pos = Poscar.from_file("CONTCAR")
                        print("copy_cmd=", copy_cmd)
                        if "ELAST" not in jobname and "LEPSILON" not in jobname:
                            # Because in ELASTIC calculations structures are deformed
                            os.system(copy_cmd)
                        # time.sleep(3)
                    except:
                        pass

            self.incar.write_file("INCAR")
            self.potcar.write_file("POTCAR")
            self.kpoints.write_file("KPOINTS")
            for i in self.copy_files:
                print("copying", i)
                shutil.copy2(i, "./")

            self.run()  # .wait()
            print("Queue 1")
            if os.path.isfile("OUTCAR"):
                try:
                    wait = Outcar(
                        "OUTCAR"
                    ).converged  # Vasprun("vasprun.xml").converged
                except:
                    pass
            print("End of the first loop", os.getcwd(), wait)

        f_energy = "na"
        enp = "na"
        contcar = str(os.getcwd()) + str("/") + str("CONTCAR")
        final_str = Poscar.from_file(contcar).atoms
        vrun = Vasprun("vasprun.xml")
        f_energy = float(vrun.final_energy)
        enp = float(f_energy) / float(final_str.num_atoms)
        natoms = final_str.num_atoms
        os.chdir("../")
        if wait == True:
            data_cal = []
            data_cal.append(
                {
                    "jobname": self.jobname,
                    "poscar": self.poscar.atoms.to_dict(),
                    "incar": self.incar.to_dict(),
                    "kpoints": self.kpoints.to_dict(),
                    "final_energy": (f_energy),
                    "contcar": final_str.to_dict(),
                }
            )
            json_file = str(self.jobname) + str(".json")
            f_json = open(json_file, "w")
            f_json.write(json.dumps(data_cal))
            f_json.close()
            print("Wrote json file", f_energy)
            return f_energy, contcar


class GenericIncars(object):
    """
    Class containing several generic Incar object for different psuedopotentials
    """

    def __init__(self, name="", incar={}, pot_type=""):
        self.name = name
        self.incar = incar
        self.pot_type = pot_type

    def optb88vdw(self):
        """
        OptB88vdW functional
        """
        data = dict(
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
            LCHARG=".FALSE.",
            LWAVE=".FALSE.",
        )
        inc = Incar(data)
        return GenericIncars(name="optb88vdw", incar=inc, pot_type="POT_GGA_PAW_PBE")

    def pbe(self):
        """
        GGA-PBE functional
        """
        data = dict(
            PREC="Accurate",
            ISMEAR=0,
            IBRION=2,
            GGA="PE",
            EDIFF="1E-7",
            NSW=1,
            NELM=400,
            ISIF=2,
            LCHARG=".FALSE.",
            LWAVE=".FALSE.",
        )
        inc = Incar(data)
        return GenericIncars(name="pbe", incar=inc, pot_type="POT_GGA_PAW_PBE")

    def lda(self):
        """
        LDA functional
        """
        data = dict(
            PREC="Accurate",
            ISMEAR=0,
            IBRION=2,
            EDIFF="1E-7",
            NSW=1,
            NELM=400,
            ISIF=2,
            LCHARG=".FALSE.",
            LWAVE=".FALSE.",
        )
        inc = Incar(data)
        return GenericIncars(name="lda", incar=inc, pot_type="POT_LDA_PAW")


"""
if __name__ == "__main__":

    p = Poscar.from_file(
        "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-RELAX-bulk@mp_541837/CONTCAR"
    )
    p = Poscar.from_file(
        "/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/CONTCAR"
    )
    # p.write_file('POSCAR')
    inc = Incar.from_file(
        "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-RELAX-bulk@mp_541837/INCAR"
    )
    inc = Incar.from_file(
        "/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/INCAR"
    )
    # inc.write_file('INCAR')
    kp = Kpoints().automatic_length_mesh(lattice_mat=p.atoms.lattice_mat)
    # kp.write_file('KPOINTS')
    # pot=Potcar(elements=p.atoms.elements)
    new_symb = list(set(p.atoms.elements))
    potcar = Potcar(elements=new_symb)
    # job = VaspJobs(poscar=p, kpoints=kp, incar=inc, jobname="testt").write_jobsub_py()
    # job = VaspJobs(poscar=p, kpoints=kp,pot_type='POT_GGA_PAW_PBE', incar=inc, jobname="testt").runjob()
    # print('optb88vdw incar',GenericIncars().optb88vdw().incar)
    # JobFactory(
    #    use_incar_dict=GenericIncars().optb88vdw().incar,
    #    pot_type=GenericIncars().optb88vdw().pot_type,
    # ).converg_encut(mat=p)
    # ).converg_kpoint(mat=p)
    # ).optimize_geometry(mat=p, encut=500, length=0)
    # ).band_structure(mat=p, encut=500, nbands=100)
    # ).loptics(mat=p, encut=500, nbands=100)
    # ).mbj_loptics(mat=p, encut=500, nbands=100)
    # ).elastic(mat=p, encut=500, nbands=100)

    JobFactory().all_optb88vdw_props(mat=p)
"""
