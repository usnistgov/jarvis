import joblib
from monty.serialization import loadfn, MontyDecoder
from jarvis.core.atoms import Atoms


def monty_to_joblib(
    read_file="/cluster/users/knc6/JARVIS-DFT1/mem-dir1L-8-22-2019.json",
    target_file="test.json",
):
    d = loadfn(read_file, cls=MontyDecoder)
    mem = []
    for i in d:
        info = {}
        jid = i["jid"]
        strt = i["data"][0]["contcar"]
        lattice = strt.lattice.matrix
        fracs = strt.frac_coords
        elements = [j.symbol for j in strt.species]
        atms = Atoms(lattice_mat=lattice, coords=fracs, elements=elements)
        phi = i["phi"]

        info["phi"] = phi
        info["atoms"] = atms
        info["jid"] = jid
        mem.append(info)
    f = open("MonolayersJOBLIB", "wb")
    joblib.dump(mem, f, protocol=2)
    f.close()


monty_to_joblib()
