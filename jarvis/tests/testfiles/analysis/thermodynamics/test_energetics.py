from jarvis.analysis.thermodynamics.energetics import (
    form_enp,
    get_twod_defect_energy,
    PhaseDiagram,
    jid_hull
)
from jarvis.db.figshare import data
import os
from jarvis.io.vasp.outputs import Vasprun
from jarvis.core.composition import Composition

tmp_xml = os.path.join(os.path.dirname(__file__), "JVASP-667_C_C_c.xml")
vrun = Vasprun(tmp_xml)

dft_3d=data('dft_3d')

def test_get_twod_defect_energy():
    Ef = get_twod_defect_energy(vrun=vrun, jid="JVASP-667", atom="C")
    print(Ef)


def test_form_enp():
    atoms = vrun.all_structures[-1]
    total_energy = vrun.final_energy
    Ef = form_enp(atoms=atoms, total_energy=total_energy)
    print(Ef)


def test_phasediag():
    system = ["Al", "O"]
    system = ["Bi", "Se"]
    system = ["Cu", "Au"]
    system = ["Al", "Ni"]
    system = ["Ga", "Al", "N"]
    system = ["Ni", "Al", "O"]
    system = ["O", "Al", "Ga", "N"]
    system = ["Ni", "Fe", "Cr"]
    system=['Mo','Se']
    x = []
    y = []
    for i in dft_3d:
        formula = i["formula"]
        comp = Composition.from_string(formula)
        atom_frac = comp.atomic_fraction
        all_elms = list(comp.to_dict())
        if (set(all_elms)).issubset(set(system)):  # and len(set(all_elms))!=1:
            x.append([i["formula"], i["formation_energy_peratom"], i["jid"]])

    pd = PhaseDiagram(
        x, verbose=True, only_plot_stable=True, only_label_stable=True
    )

    pd.plot()
def test_jid_hull():
   x=jid_hull(jid='JVASP-114723',dataset=dft_3d)
# test_get_twod_defect_energy()
# test_form_enp()
