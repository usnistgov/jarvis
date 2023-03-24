"""Module to converge K-points and Cutoff."""
from jarvis.core.atoms import Atoms
from jarvis.core.kpoints import Kpoints3D as Kpoints
from jarvis.tasks.qe.qe import QEjob
import time

sleep = 5


def converg_kpoints(
    atoms=None,
    length=0,
    encut=40,
    ecutrho=250,
    tol=0.001,
    increment=5,
    qe_cmd="/cluster/deb9/bin/mpirun -n 16 /cluster/bin/pw.x",
    psp_dir=None,
    url=None,
    psp_temp_name=None,
):
    """Converge k-points for a material."""
    scf_init = {
        "control": {
            "calculation": "'scf'",
            "restart_mode": "'from_scratch'",
            "prefix": "'KPOINTS'",
            "outdir": "'./'",
            "tstress": ".true.",
            "tprnfor": ".true.",
            "disk_io": "'nowf'",
            "pseudo_dir": None,
            "verbosity": "'high'",
            "nstep": 100,
            "etot_conv_thr": "1.0d-5",
        },
        "system": {
            "ibrav": 0,
            "degauss": 0.01,
            "nat": None,
            "ntyp": None,
            "ecutwfc": encut,
            "ecutrho": ecutrho,
            "occupations": "'smearing'",
            "smearing": "'mp'",
        },
        "electrons": {
            "diagonalization": "'david'",
            "mixing_mode": "'plain'",
            "mixing_beta": 0.7,
            "conv_thr": "1d-9",
        },
    }

    en1 = -10000
    convg_kp1 = False
    convg_kp2 = False
    length1 = length
    kp_list = []
    while not convg_kp2 and not convg_kp1:
        length1 = length1 + increment
        print("Incrementing length", length1)
        kpoints = Kpoints().automatic_length_mesh(
            lattice_mat=atoms.lattice_mat, length=length
        )

        prefix = "KPOINTS-" + str(length1)
        scf_init["control"]["prefix"] = str('"') + prefix + str('"')
        kpoints = Kpoints().automatic_length_mesh(
            lattice_mat=atoms.lattice_mat, length=length
        )
        mesh = kpoints.kpts[0]
        if mesh not in kp_list:
            kp_list.append(mesh)
            scf_init["system"]["ntyp"] = ""
            qejob_scf_init = QEjob(
                atoms=atoms,
                psp_dir=psp_dir,
                input_params=scf_init,
                output_file="scf_init.out",
                qe_cmd=qe_cmd,
                jobname=prefix,
                kpoints=kpoints,
                psp_temp_name=psp_temp_name,
                url=url,
                input_file=prefix + "_ascf_init.in",
            )
            info_scf = qejob_scf_init.runjob()
            print("Energy", info_scf["total_energy"])
            en2 = float(info_scf["total_energy"])
            time.sleep(sleep)
            while abs(en2 - en1) > tol:
                en1 = en2
                print("Incrementing length", length1)
                while mesh in kp_list:
                    length1 = length1 + increment
                    # Assuming you are not super unlucky
                    # kpoints = Auto_Kpoints(mat=mat, length=length1)
                    kpoints = Kpoints().automatic_length_mesh(
                        lattice_mat=atoms.lattice_mat, length=length1
                    )  # Auto_Kpoints(mat=mat, length=length)
                    mesh = kpoints.kpts[0]

                kpoints = Kpoints().automatic_length_mesh(
                    lattice_mat=atoms.lattice_mat, length=length1
                )  # Auto_Kpoints(mat=mat, length=length)
                mesh = kpoints.kpts[0]
                prefix = "KPOINTS-" + str(length1)
                scf_init["control"]["prefix"] = str('"') + prefix + str('"')
                if mesh not in kp_list:
                    kp_list.append(mesh)

                    scf_init["system"]["ntyp"] = ""
                    qejob_scf_init = QEjob(
                        atoms=atoms,
                        input_params=scf_init,
                        psp_dir=psp_dir,
                        output_file="scf_init.out",
                        qe_cmd=qe_cmd,
                        jobname=prefix,
                        kpoints=kpoints,
                        psp_temp_name=psp_temp_name,
                        url=url,
                        input_file=prefix + "_ascf_init.in",
                    )

                    info_scf = qejob_scf_init.runjob()
                    en2 = float(info_scf["total_energy"])
                    print("Energy", info_scf["total_energy"])
                    time.sleep(sleep)
                else:

                    length1 = length1 + increment
                    # Assuming you are not super unlucky
                    # kpoints = Auto_Kpoints(mat=mat, length=length1)
                    kpoints = Kpoints().automatic_length_mesh(
                        lattice_mat=atoms.lattice_mat, length=length1
                    )  # Auto_Kpoints(mat=mat, length=length)
                    mesh = kpoints.kpts[0]
                    kp_list.append(mesh)

            convg_kp1 = True

            # Some extra points to check
            print("Some extra points to check for KPOINTS")
            length3 = length1 + increment
            prefix = "KPOINTS-" + str(length3)
            scf_init["control"]["prefix"] = str('"') + prefix + str('"')

            kpoints = Kpoints().automatic_length_mesh(
                lattice_mat=atoms.lattice_mat, length=length3
            )
            mesh = kpoints.kpts[0]
            kp_list.append(mesh)
            scf_init["system"]["ntyp"] = ""
            qejob_scf_init = QEjob(
                atoms=atoms,
                psp_dir=psp_dir,
                input_params=scf_init,
                output_file="scf_init.out",
                qe_cmd=qe_cmd,
                psp_temp_name=psp_temp_name,
                url=url,
                jobname=prefix,
                kpoints=kpoints,
                input_file=prefix + "_ascf_init.in",
            )
            info_scf = qejob_scf_init.runjob()
            en3 = float(info_scf["total_energy"])
            print("Energy", info_scf["total_energy"])
            time.sleep(sleep)

            length4 = length3 + increment
            prefix = "KPOINTS-" + str(length4)
            scf_init["control"]["prefix"] = str('"') + prefix + str('"')

            kpoints = Kpoints().automatic_length_mesh(
                lattice_mat=atoms.lattice_mat, length=length4
            )
            mesh = kpoints.kpts[0]
            kp_list.append(mesh)
            scf_init["system"]["ntyp"] = ""
            qejob_scf_init = QEjob(
                atoms=atoms,
                psp_dir=psp_dir,
                input_params=scf_init,
                output_file="scf_init.out",
                qe_cmd=qe_cmd,
                psp_temp_name=psp_temp_name,
                url=url,
                jobname=prefix,
                kpoints=kpoints,
                input_file=prefix + "_ascf_init.in",
            )
            info_scf = qejob_scf_init.runjob()
            en4 = info_scf["total_energy"]
            print("Energy", info_scf["total_energy"])
            time.sleep(sleep)

            length5 = length4 + increment
            prefix = "KPOINTS-" + str(length5)
            scf_init["control"]["prefix"] = str('"') + prefix + str('"')
            kpoints = Kpoints().automatic_length_mesh(
                lattice_mat=atoms.lattice_mat, length=length5
            )
            mesh = kpoints.kpts[0]
            kp_list.append(mesh)
            scf_init["system"]["ntyp"] = ""
            qejob_scf_init = QEjob(
                atoms=atoms,
                psp_dir=psp_dir,
                input_params=scf_init,
                output_file="scf_init.out",
                qe_cmd=qe_cmd,
                jobname=prefix,
                kpoints=kpoints,
                psp_temp_name=psp_temp_name,
                url=url,
                input_file=prefix + "_ascf_init.in",
            )
            info_scf = qejob_scf_init.runjob()
            en5 = info_scf["total_energy"]
            print("Energy", info_scf["total_energy"])
            time.sleep(sleep)

            length6 = length5 + increment
            prefix = "KPOINTS-" + str(length6)
            scf_init["control"]["prefix"] = str('"') + prefix + str('"')
            kpoints = Kpoints().automatic_length_mesh(
                lattice_mat=atoms.lattice_mat, length=length6
            )
            mesh = kpoints.kpts[0]
            kp_list.append(mesh)
            scf_init["system"]["ntyp"] = ""
            qejob_scf_init = QEjob(
                atoms=atoms,
                input_params=scf_init,
                psp_dir=psp_dir,
                output_file="scf_init.out",
                qe_cmd=qe_cmd,
                jobname=prefix,
                kpoints=kpoints,
                psp_temp_name=psp_temp_name,
                url=url,
                input_file=prefix + "_ascf_init.in",
            )
            info_scf = qejob_scf_init.runjob()
            en6 = info_scf["total_energy"]
            print("Energy", info_scf["total_energy"])
            time.sleep(sleep)

            length7 = length6 + increment
            prefix = "KPOINTS-" + str(length7)
            scf_init["control"]["prefix"] = str('"') + prefix + str('"')
            kpoints = Kpoints().automatic_length_mesh(
                lattice_mat=atoms.lattice_mat, length=length7
            )
            mesh = kpoints.kpts[0]
            kp_list.append(mesh)
            scf_init["system"]["ntyp"] = ""
            qejob_scf_init = QEjob(
                atoms=atoms,
                input_params=scf_init,
                output_file="scf_init.out",
                psp_dir=psp_dir,
                qe_cmd=qe_cmd,
                jobname=prefix,
                kpoints=kpoints,
                psp_temp_name=psp_temp_name,
                url=url,
                input_file=prefix + "_ascf_init.in",
            )
            info_scf = qejob_scf_init.runjob()
            en7 = info_scf["total_energy"]
            print("Energy", info_scf["total_energy"])
            time.sleep(sleep)

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
                print(
                    "KPOINTS convergence achieved for ",
                    length1,
                )
                convg_kp2 = True

            return length1


def converg_cutoff(
    atoms=None,
    length=10,
    encut=40,
    ecutrho=250,
    tol=0.001,
    increment=5,
    psp_dir=None,
    url=None,
    psp_temp_name=None,
    qe_cmd="/cluster/deb9/bin/mpirun -n 16 /cluster/bin/pw.x",
):
    """Converge cutoff for a material."""
    scf_init = {
        "control": {
            "calculation": "'scf'",
            "restart_mode": "'from_scratch'",
            "prefix": "'ENCUT'",
            "outdir": "'./'",
            "tstress": ".true.",
            "tprnfor": ".true.",
            "disk_io": "'nowf'",
            "pseudo_dir": None,
            "verbosity": "'high'",
            "nstep": 100,
            "etot_conv_thr": "1.0d-5",
        },
        "system": {
            "ibrav": 0,
            "degauss": 0.01,
            "nat": None,
            "ntyp": None,
            "ecutwfc": encut,
            "ecutrho": ecutrho,
            "occupations": "'smearing'",
            "smearing": "'mp'",
        },
        "electrons": {
            "diagonalization": "'david'",
            "mixing_mode": "'plain'",
            "mixing_beta": 0.7,
            "conv_thr": "1d-9",
        },
    }

    en1 = -10000
    convg_encut1 = False
    convg_encut2 = False
    encut1 = encut
    kpoints = Kpoints().automatic_length_mesh(
        lattice_mat=atoms.lattice_mat, length=length
    )
    while not convg_encut2 and not convg_encut1:

        prefix = "ENCUT-" + str(encut1)
        scf_init["control"]["prefix"] = str('"') + prefix + str('"')
        scf_init["system"]["ntyp"] = ""
        scf_init["system"]["ecutwfc"] = encut
        qejob_scf_init = QEjob(
            atoms=atoms,
            psp_dir=psp_dir,
            input_params=scf_init,
            output_file="scf_init.out",
            qe_cmd=qe_cmd,
            jobname=prefix,
            kpoints=kpoints,
            url=url,
            psp_temp_name=psp_temp_name,
            input_file=prefix + "_ascf_init.in",
        )
        info_scf = qejob_scf_init.runjob()
        print("Energy", info_scf["total_energy"])
        en2 = float(info_scf["total_energy"])
        time.sleep(sleep)
        while abs(en2 - en1) > tol:
            en1 = en2
            encut1 = encut1 + increment
            print("Incrementing cutoff", encut1)
            prefix = "ENCUT-" + str(encut1)
            scf_init["control"]["prefix"] = str('"') + prefix + str('"')
            scf_init["system"]["ntyp"] = ""
            scf_init["system"]["ecutwfc"] = encut1
            qejob_scf_init = QEjob(
                atoms=atoms,
                input_params=scf_init,
                psp_dir=psp_dir,
                output_file="scf_init.out",
                qe_cmd=qe_cmd,
                jobname=prefix,
                kpoints=kpoints,
                url=url,
                psp_temp_name=psp_temp_name,
                input_file=prefix + "_ascf_init.in",
            )

            info_scf = qejob_scf_init.runjob()
            en2 = float(info_scf["total_energy"])
            print("Energy", info_scf["total_energy"])
            time.sleep(sleep)

        convg_encut1 = True

        # Some extra points to check
        print("Some extra points to check for ENCUT")
        encut3 = encut1 + increment
        prefix = "ENCUT-" + str(encut3)
        scf_init["control"]["prefix"] = str('"') + prefix + str('"')
        scf_init["system"]["ntyp"] = ""
        scf_init["system"]["ecutwfc"] = encut3

        qejob_scf_init = QEjob(
            atoms=atoms,
            input_params=scf_init,
            output_file="scf_init.out",
            qe_cmd=qe_cmd,
            psp_dir=psp_dir,
            jobname=prefix,
            kpoints=kpoints,
            url=url,
            psp_temp_name=psp_temp_name,
            input_file=prefix + "_ascf_init.in",
        )
        info_scf = qejob_scf_init.runjob()
        en3 = float(info_scf["total_energy"])
        print("Energy", info_scf["total_energy"])
        time.sleep(sleep)

        encut4 = encut3 + increment
        prefix = "ENCUT-" + str(encut4)
        scf_init["control"]["prefix"] = str('"') + prefix + str('"')
        scf_init["system"]["ntyp"] = ""
        scf_init["system"]["ecutwfc"] = encut4
        qejob_scf_init = QEjob(
            atoms=atoms,
            input_params=scf_init,
            psp_dir=psp_dir,
            output_file="scf_init.out",
            qe_cmd=qe_cmd,
            jobname=prefix,
            kpoints=kpoints,
            url=url,
            psp_temp_name=psp_temp_name,
            input_file=prefix + "_ascf_init.in",
        )
        info_scf = qejob_scf_init.runjob()
        en4 = info_scf["total_energy"]
        print("Energy", info_scf["total_energy"])
        time.sleep(sleep)

        encut5 = encut4 + increment
        prefix = "ENCUT-" + str(encut5)
        scf_init["control"]["prefix"] = str('"') + prefix + str('"')
        scf_init["system"]["ntyp"] = ""
        scf_init["system"]["ecutwfc"] = encut5
        qejob_scf_init = QEjob(
            atoms=atoms,
            input_params=scf_init,
            output_file="scf_init.out",
            qe_cmd=qe_cmd,
            psp_dir=psp_dir,
            jobname=prefix,
            kpoints=kpoints,
            url=url,
            psp_temp_name=psp_temp_name,
            input_file=prefix + "_ascf_init.in",
        )
        info_scf = qejob_scf_init.runjob()
        en5 = info_scf["total_energy"]
        print("Energy", info_scf["total_energy"])
        time.sleep(sleep)

        encut6 = encut5 + increment
        prefix = "ENCUT-" + str(encut6)
        scf_init["control"]["prefix"] = str('"') + prefix + str('"')
        scf_init["system"]["ntyp"] = ""
        scf_init["system"]["ecutwfc"] = encut6
        qejob_scf_init = QEjob(
            atoms=atoms,
            input_params=scf_init,
            output_file="scf_init.out",
            psp_dir=psp_dir,
            qe_cmd=qe_cmd,
            jobname=prefix,
            kpoints=kpoints,
            url=url,
            psp_temp_name=psp_temp_name,
            input_file=prefix + "_ascf_init.in",
        )
        info_scf = qejob_scf_init.runjob()
        en6 = info_scf["total_energy"]
        print("Energy", info_scf["total_energy"])
        time.sleep(sleep)

        encut7 = encut6 + increment
        prefix = "ENCUT-" + str(encut7)
        scf_init["control"]["prefix"] = str('"') + prefix + str('"')
        scf_init["system"]["ntyp"] = ""
        scf_init["system"]["ecutwfc"] = encut7
        qejob_scf_init = QEjob(
            atoms=atoms,
            input_params=scf_init,
            output_file="scf_init.out",
            qe_cmd=qe_cmd,
            jobname=prefix,
            psp_dir=psp_dir,
            kpoints=kpoints,
            url=url,
            psp_temp_name=psp_temp_name,
            input_file=prefix + "_ascf_init.in",
        )
        info_scf = qejob_scf_init.runjob()
        en7 = info_scf["total_energy"]
        print("Energy", info_scf["total_energy"])
        time.sleep(sleep)

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
            print("ENCUT convergence achieved for ", encut)
            convg_encut2 = True
        return encut


if __name__ == "__main__":
    from jarvis.db.figshare import get_jid_data

    atoms = Atoms.from_dict(
        get_jid_data(jid="JVASP-1002", dataset="dft_3d")["atoms"]
    )
    # converg_kpoints(atoms=atoms)
    converg_cutoff(atoms=atoms)
