import os
from jarvis.io.vasp.outputs import Oszicar
from jarvis.io.vasp.inputs import Poscar, Incar
from jarvis.core.kpoints import Kpoints3D as Kpoints


class VaspJobErrorHandlers(object):
    def __init__(self, run_directory, jobid="", all_errors=[]):
        self.run_directory = run_directory
        self.jobid = jobid
        self.all_errors = all_errors

    def check_errors(self, logfile="vasp.out", timeout=18000):
        os.chdir(self.run_directory)

        errors = []
        # print ("going here 12")
        # try:
        #  run=Vasprun("vasprun.xml")
        #  if run.converged_electronic == False:
        #      errors.append("unconverged_electronic")
        #  if run.converged_ionic == False:
        #      errors.append("unconverged")
        # except:
        #    pass
        # f=open("OUTCAR",'r')
        try:
            oszicar = Oszicar("OSZICAR")
            ionic_steps = len(oszicar.ionic_steps)
            with open("OUTCAR") as f:
                for line in f:

                    if "NSW" in line:
                        try:
                            # nbands = int(d[-1].strip())
                            nsw = int(
                                str(line.split("=")[1]).split("number of steps")[0]
                            )
                            if nsw > 1 and ionic_steps == nsw:
                                errors.append("unconverged")
                        except:
                            pass

                    if "NELM" in line:
                        try:
                            # nbands = int(d[-1].strip())
                            # nsw= int(str(line.split("=")[1]).split("number of steps")[0])
                            nelm = int(str(line.split("=")[1]).split(";")[0])
                            electronic_steps = int(oszicar.electronic_steps[-1][1])
                            if electronic_steps < nelm:
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
                if "WARNING: Sub-Space-Matrix is not hermitian in" in line:
                    err = "subspacematrix"
                    if err not in errors:
                        errors.append(err)
                if "Tetrahedron method fails for NKPT<" in line:
                    err = "tet"
                    if err not in errors:
                        errors.append(err)
                if "Fatal error detecting k-mesh" in line:
                    err = "tet"
                    if err not in errors:
                        errors.append(err)
                if "Fatal error: unable to match k-point" in line:
                    err = "tet"
                    if err not in errors:
                        errors.append(err)
                if "Routine TETIRR needs special values" in line:
                    err = "tet"
                    if err not in errors:
                        errors.append(err)
                if "Tetrahedron method fails for NKPT<" in line:
                    err = "tet"
                    if err not in errors:
                        errors.append(err)
                if "inverse of rotation matrix was not found (increase" in line:
                    err = "inv_rot_ma"
                    if err not in errors:
                        errors.append(err)
                if "SYMPREC" in line:
                    err = "inv_rot_ma"
                    if err not in errors:
                        errors.append(err)
                if "Routine TETIRR needs special values" in line:
                    err = "tetirr"
                    if err not in errors:
                        errors.append(err)
                if "Could not get correct shift" in line:
                    err = "incorrect_shift"
                    if err not in errors:
                        errors.append(err)
                if "REAL_OPTLAY: internal error" in line:
                    err = "real_optlay"
                    if err not in errors:
                        errors.append(err)
                if "REAL_OPT: internal ERROR" in line:
                    err = "real_optlay"
                    if err not in errors:
                        errors.append(err)
                if "ERROR RSPHER" in line:
                    err = "rspher"
                    if err not in errors:
                        errors.append(err)
                if "DENTET" in line:
                    err = "dentet"
                    if err not in errors:
                        errors.append(err)
                if "TOO FEW BAND" in line:
                    err = "too_few_bands"
                    if err not in errors:
                        errors.append(err)
                if "ERROR: the triple product of the basis vectors" in line:
                    err = "triple_product"
                    if err not in errors:
                        errors.append(err)
                if "Found some non-integer element in rotation matrix" in line:
                    err = "rot_matrix"
                    if err not in errors:
                        errors.append(err)
                if "BRIONS problems: POTIM should be increased" in line:
                    err = "brions"
                    if err not in errors:
                        errors.append(err)
                if "internal error in subroutine PRICEL" in line:
                    err = "pricel"
                    if err not in errors:
                        errors.append(err)
                if "LAPACK: Routine ZPOTRF failed" in line:
                    err = "zpotrf"
                    if err not in errors:
                        errors.append(err)
                if "One of the lattice vectors is very long (>50 A), but AMIN" in line:
                    err = "amin"
                    if err not in errors:
                        errors.append(err)
                if "ZBRENT: fatal internal in" in line:
                    err = "zbrent"
                    if err not in errors:
                        errors.append(err)
                if "ZBRENT: fatal error in bracketing" in line:
                    err = "zrbent"
                    if err not in errors:
                        errors.append(err)

                if "ERROR in subspace rotation PSSYEVX" in line:
                    err = "pssyevx"
                    if err not in errors:
                        errors.append(err)
                if "WARNING in EDDRMM: call to ZHEGV failed" in line:
                    err = "eddrmm"
                    if err not in errors:
                        errors.append(err)
                if "Error EDDDAV: Call to ZHEGV failed" in line:
                    err = "edddav"
                    if err not in errors:
                        errors.append(err)
                if "Your FFT grids (NGX,NGY,NGZ) are not sufficient" in line:
                    err = "aliasing_incar"
                    if err not in errors:
                        errors.append(err)

        # UNCONVERGED
        st = os.stat(logfile)
        if time.time() - st.st_mtime > timeout:
            errors.append("Frozenjob")

        all_errors = list(set(errors))
        self.all_errors = all_errors
        return all_errors

    def correct_errors(self):
        msg = self.all_errors
        mat = Poscar.from_file("POSCAR")
        incar = Incar.from_file("INCAR")
        kpoints = Kpoints.from_file("KPOINTS")
        while len(msg) > 0:

            if "zpotrf" in msg:
                print("zpotrf error")

                try:
                    oszicar = Oszicar("OSZICAR")
                    nsteps = len(oszicar.ionic_steps)
                except:
                    nsteps = 0

                if nsteps >= 1:
                    potim = float(incar.get("POTIM", 0.5)) / 2.0
                    incar.update({"ISYM": 0, "POTIM": potim})
                    # incar.write_file("INCAR")
                else:
                    s = mat.atoms
                    comm = mat.comment
                    s.apply_strain(0.05)
                    pos = Poscar(s)
                    pos.comment = comm + str("corrected")
                    # line = str(pos.comment) + str(" ") + str(s) + "\n"
                    # correc.write(line)
                    # pos.write_file("POSCAR")
                    try:
                        cmd = str("rm -rf WAVECAR CHGCAR")
                        os.system(cmd)
                    except:
                        pass

            # if "brions" in msg:
            #     potim = float(incar.get("POTIM", 0.5)) +0.1
            #     incar.update({"POTIM":potim})

            if "zbrent" in msg:
                incar.update({"IBRION": 1})
            if "too_few_bands" in msg:
                with open("OUTCAR") as f:
                    for line in f:
                        if "NBANDS" in line:
                            try:
                                d = line.split("=")
                                nbands = int(d[-1].strip())
                                incar.update({"NBANDS": int(1.1 * nbands)})
                                break
                            except (IndexError, ValueError):
                                pass

            if "pssyevx" in msg:
                incar.update({"ALGO": "Normal"})
            if "edddav" in msg:
                incar.update({"LPLANE": "False"})
                incar.update({"ALGO": "Damped"})
                try:
                    cmd = str("rm -rf CHGCAR")
                    os.system(cmd)
                except:
                    pass
            if "pricel" in msg:
                incar.update({"ALGO": "Normal"})
            if "subspacematrix" in msg:
                incar.update({"NPAR": 1})
            if "rspher" in msg:
                incar.update({"NPAR": 1})
            if "real_optlay" in msg:
                incar.update({"NPAR": 1})
            if "subspacematrix" in msg:
                incar.update({"NPAR": 1})

            if "rot_matrix" in msg:
                incar.update({"ISYM": 0})
            if "amin" in msg:
                incar.update({"AMIN": 0.01})
            if "FrozenJob" in msg:
                try:
                    cmd = str("qdel ") + str(self.job_id)
                    os.system(cmd)
                except:
                    pass

                try:
                    cmd = str("scancel ") + str(self.job_id)
                    os.system(cmd)
                except:
                    pass
                # incar.update({"EDIFF":1e-8,"NSW":1,"IBRION":1})
            # Unconverged
            if "unconverged_electronic" in msg:
                nelm = int(incar.get("NELM", 60))
                incar.update({"NELM": int(nelm) + 100})
            incar.write_file("INCAR")
            kpoints.write_file("KPOINTS")
            potcar.write_file("POTCAR")
            mat.write_file("POSCAR")
            if "unconverged" in msg and os.path.isfile("OSZICAR"):
                oszicar = Oszicar("OSZICAR")
                nsteps = len(oszicar.ionic_steps)

                with open("OUTCAR") as f:
                    for line in f:

                        if "NSW" in line:
                            try:
                                nsw = int(
                                    str(line.split("=")[1]).split("number of steps")[0]
                                )
                                incar.update({"EDIFF": 1e-8, "NSW": nsw})
                                shutil.copy2("POSCAR", "POSCAR.orig")
                                # shutil.copy2('CONTCAR','POSCAR')
                            except:
                                pass
                if nsw > 1 and nsteps == nsw:
                    shutil.copy2("CONTCAR", "POSCAR")
