from jarvis.io.vasp.inputs import Poscar
import os
import json

# OptB88vdW energy per atoms for elements
unary_json_file = str(os.path.join(os.path.dirname(__file__), "unary.json"))
unary_json = open(unary_json_file, "r")
unary_data = json.load(unary_json)
unary_json.close()


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def unary_energy(el="Na"):
    en = "na"
    for i, j in unary_data.items():
        if str(i) == str(el):
            en = j["energy"]
    return en


def form_enp(atoms=None, total_energy=None):
    """
    Calculate formation energy given the total energy and the atoms object
    Currently for OptB88vdW functional based chemical potential implemented
    but can be generalized by changing unary_energy
    """
    dd = atoms.composition.to_dict()
    ref = 0.0
    for k, v in dd.items():
        e1 = unary_energy(k)
        if e1 == "na":
            ref = "na"
            ValueError("Element reference does not exist", el)

        else:
            ref = ref + float(v) * float(e1)
    if isfloat(ref):
        form_en = float(total_energy) - float(ref)
        form_en = round(float(form_en) / float(atoms.num_atoms), 5)
    return form_en


# https://wiki.fysik.dtu.dk/ase/_modules/ase/phasediagram.html#PhaseDiagram
# class EnergyConvexHull(object):
#     def __init__(entries=None):
#        # [array of [composition_as_dict,total_energy]]
#        self.entries = entries

#     def chull(self):
#        for i in self.entries:

"""
if __name__ == "__main__":
    p = Poscar.from_file(
        "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-RELAX-bulk@mp_541837/POSCAR"
    ).atoms
    total_energy = -9.974648
    fen = form_enp(atoms=p, total_energy=total_energy)
    print("fen", fen)
    x = [
        {"comp": {"Al": 2, "O": 1}, "energy": -0.15, "id": "xyz"},
        {"comp": {"Al": 2, "O": 3}, "energy": -16.5, "id": "xyz"},
        {"comp": {"Al": 1, "O": 1}, "energy": -4.86, "id": "xyz"},
    ]

    # entries = [[{'Al':2,'O':1},-0.15],[{'Al':2,'O':3},-16.065],[{'Al':1,'O':1},-2.8],[{'Al':1,'O':3},-4.865]]
    print(x)
"""
