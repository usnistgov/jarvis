from __future__ import absolute_import
import os
from jarvis.core.atoms import Atoms
from jarvis.io.phonopy.outputs import (
    bandstructure_plot,
    total_dos,
    read_fc,
    get_phonon_tb,
)
from jarvis.io.phonopy.inputs import PhonopyInputs
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.io.wannier.outputs import WannierHam
import requests, tempfile, zipfile, io
from jarvis.io.vasp.outputs import Vasprun

# from phonopy.interface.vasp import read_vasp
from jarvis.analysis.phonon.force_constants import qpoint, read_fc
from jarvis.io.phonopy.outputs import get_phonon_tb
import matplotlib.pyplot as plt

plt.switch_backend("agg")
from jarvis.db.figshare import data

fls = data("raw_files")
fc_file = os.path.join(os.path.dirname(__file__), "FORCE_CONSTANTS")
pos = os.path.join(os.path.dirname(__file__), "POSCAR")
wtb = os.path.join(os.path.dirname(__file__), "phonopyTB_hr.dat")
totdos = os.path.join(os.path.dirname(__file__), "total_dos.dat")
band = os.path.join(os.path.dirname(__file__), "band.yaml")


def test_fc():
    fc = read_fc(fc_file)
    assert (fc[0][0][0][0]) == 12.960974735
    qp = qpoint(force_constant=fc, qpt=[0, 0, 0])
    assert (qp[0][0][0][0]) == 12.960974735


def test_wann():
    a = Atoms.from_poscar(pos)
    fc = read_fc(fc_file)
    get_phonon_tb(fc=fc, atoms=a, out_file=wtb)
    cvn = Spacegroup3D(a).conventional_standard_structure
    w = WannierHam(wtb)
    w.get_bandstructure_plot(atoms=cvn, yrange=[0, 550])
    cmd = "rm phonopyTB_hr.dat bs.png"
    os.system(cmd)


def test_outputs():
    # JVASP-1002
    x, y = total_dos(tot_dos=totdos)
    a, b, c, d = bandstructure_plot(band_yaml=band)


def test_download(jid="JVASP-1002"):
    for i in fls["FD-ELAST"]:
        if isinstance(i, dict):
            if i["name"].split(".zip")[0] == jid:
                print(i)

                r = requests.get(i["download_url"])
                z = zipfile.ZipFile(io.BytesIO(r.content))
                vrun_path = z.read("vasprun.xml").decode("utf-8")
                fd, path = tempfile.mkstemp()
                with os.fdopen(fd, "w") as tmp:
                    tmp.write(vrun_path)
                vrun = Vasprun(path)
                fc = vrun.phonon_data()["force_constants"]
                atoms = vrun.all_structures[0]
                print(atoms)
                # atoms = Atoms.from_poscar(pos)
                print(atoms)
                fd, path = tempfile.mkstemp()
                get_phonon_tb(fc=fc, atoms=atoms, out_file=path)
                cvn = Spacegroup3D(atoms).conventional_standard_structure
                w = WannierHam(path)
                # print ('atoms',atoms)
                w.get_bandstructure_plot(atoms=atoms, yrange=[0, 550])


test_download()
