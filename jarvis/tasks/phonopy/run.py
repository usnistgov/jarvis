"""Module to run phonopy. Assumes phonopy in PATH."""
from jarvis.io.vasp.outputs import Vasprun
import os
from jarvis.io.phonopy.inputs import PhonopyInputs


def run_phonopy(path="MAIN-ELAST"):
    """Run phonopy."""
    os.chdir(path)
    try:
        vrun = Vasprun(os.path.join(path, "vasprun.xml"))
        atoms = vrun.all_Structures[-1]
        PhonopyInputs(atoms=atoms).generate_all_files()
        cmd = "phonopy -p meshdos.conf"
        os.system(cmd)
        cmd = "phonopy -t meshdos.conf"
        os.system(cmd)
        cmd = "phonopy -p bandd.conf"
        os.system(cmd)
        cmd = "bandplot -o PBAND.png"
        os.system(cmd)
    except Exception:
        print("Cannot run phonopy. Check phonopy PATH and vasprun.")
        pass
