from jarvis.db.static.explore_db import get_3d_dataset
import os
from pymatgen.symmetry.bandstructure import HighSymmKpath

data = get_3d_dataset()


def give_cif_for_jid(jid=""):
    strt = "na"
    for i in data:
        if i["jid"] == jid:
            strt = i["final_str"]
            break
    return strt


def prepare_wien_input(jid="JVASP-1002"):
    s = give_cif_for_jid(jid=jid)
    filename1 = str(jid) + str(".cif")
    fold = jid

    ## if not os.path.exists(fold):
    ##os.makedirs(fold)
    ##os.chdir(fold)
    s.to(fmt="cif", filename=filename1)
    cmd = str("cif2struct") + str(" ") + str(filename1)
    os.system(cmd)
    cmd = str("init_lapw -b -red 0 -vxc 13 -ecut -7.0  -numk 100 -sp >init_w2k_out")
    os.system(cmd)
    cmd = str("runsp_lapw -cc 0.0001 -ec 0.0001 -p -i 500  >run_w2k_out")
    os.system(cmd)
    ##os.chdir('../')


def get_kpoints(strt=""):
    kpath = HighSymmKpath(strt)
    kp = kpath.kpath["kpoints"]
    uniqe_lbl = []
    uniqe_k = []
    for i, j in kp.items():
        uniqe_lbl.append(j)
        uniqe_k.append(i)
    return uniqe_lbl, uniqe_k


"""
from pymatgen.core.structure import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
s=Structure.from_file('JVASP-1002.cif')
x,y=get_kpoints(s) 
for i,j in zip(x,y):
  print (i,j)
"""
