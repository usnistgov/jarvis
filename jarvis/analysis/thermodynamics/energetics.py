"""Get formation energy."""

import os
from jarvis.db.figshare import data
from jarvis.core.atoms import Atoms
from jarvis.db.jsonutils import loadjson


def get_optb88vdw_energy():
    """Get OptB88vdW energy per atoms for elements."""
    return loadjson(os.path.join(os.path.dirname(__file__), "unary.json"))


def get_unary_qe_tb_energy():
    """Get elemental chemical potential for GBRV tight-binding project."""
    return loadjson(
        os.path.join(os.path.dirname(__file__), "unary_qe_tb.json")
    )


def isfloat(value):
    """Check if a number is float.

    TODO: replace with isinstance.
    """
    try:
        float(value)
        return True
    except ValueError:
        return False


def unary_energy(el="Na", chem_pots=None):
    """Provide energy per atoms of an element."""
    if chem_pots is None:
        chem_pots = get_optb88vdw_energy()
    en = "na"
    for i, j in chem_pots.items():
        if str(i) == str(el):
            en = j["energy"]
    return en


def form_enp(atoms=None, total_energy=None, chem_pots=None):
    """
    Calculate formation energy given the total energy and the atoms object.

    Currently for OptB88vdW functional based chemical potential implemented
    but can be generalized by changing unary_energy.
    """
    dd = atoms.composition.to_dict()
    # print ('dd',dd)
    ref = 0.0
    for k, v in dd.items():
        e1 = unary_energy(el=k, chem_pots=chem_pots)
        # print (k,v,e1,total_energy)
        if e1 == "na":
            ref = "na"
            ValueError("Element reference does not exist", e1)

        else:
            ref = ref + float(v) * float(e1)
    if isfloat(ref):
        form_en = float(total_energy) - float(ref)
        form_en = round(float(form_en) / float(atoms.num_atoms), 5)
    return form_en


def get_twod_defect_energy(vrun="", jid="", atom=""):
    """Get mono 2D defect formation energy with OptB88vdW data."""
    dft2d = data("dft_2d")

    def get_enp_jid(jid=""):
        for i in dft2d:
            if i["jid"] == jid:
                return (
                    i["optb88vdw_total_energy"]
                    / Atoms.from_dict(i["atoms"]).num_atoms
                )

        # dir='JVASP-667_C_C_c'
        # tmp=dir.split('_')
        # jid=tmp[0]
        # atom=tmp[2]

    strt = vrun.all_structures[-1]
    natoms = strt.num_atoms
    fin_en = vrun.final_energy
    chem_pot = unary_energy(atom)
    bulk_en_pa = get_enp_jid(jid)
    Ef = fin_en - (natoms + 1) * bulk_en_pa + chem_pot
    return Ef


# https://wiki.fysik.dtu.dk/ase/_modules/ase/phasediagram.html#PhaseDiagram
# class EnergyConvexHull(object):
#     def __init__(entries=None):
#        # [array of [composition_as_dict,total_energy]]
#        self.entries = entries

#     def chull(self):
#        for i in self.entries:

"""
if __name__ == "__main__":
    from jarvis.io.vasp.inputs import Poscar
    p = Poscar.from_file(
        "../..//POSCAR"
    ).atoms
    total_energy = -9.974648
    fen = form_enp(atoms=p, total_energy=total_energy)
    print("fen", fen)
    x = [
        {"comp": {"Al": 2, "O": 1}, "energy": -0.15, "id": "xyz"},
        {"comp": {"Al": 2, "O": 3}, "energy": -16.5, "id": "xyz"},
        {"comp": {"Al": 1, "O": 1}, "energy": -4.86, "id": "xyz"},
    ]

    # entries = [[{'Al':2,'O':1},-0.15],[{'Al':2,'O':3},
     -16.065],[{'Al':1,'O':1},-2.8],[{'Al':1,'O':3},-4.865]]
    print(x)
"""
