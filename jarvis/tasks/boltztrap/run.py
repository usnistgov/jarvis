"""Module to run BoltzTrap. Assumes x_trans in PATH."""
import os
from jarvis.io.boltztrap.inputs import WriteInputs


def run_boltztrap(
    path="MAIN-RELAX", cmd="~/anaconda2/bin/x_trans BoltzTraP -so"
):
    """Run boltztrap program with vasprun.xml."""
    os.chdir(path)
    try:
        vrun = os.path.join(path, "vasprun.xml")
        inp = WriteInputs(vasprun_path=vrun)
        inp.write_energy()
        inp.write_struct()
        inp.write_intrans()
        os.system(cmd)
    except Exception:
        print("Error running boltztrap, check x_trans and vasprun.")
        pass
    os.chdir(path)
