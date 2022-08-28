from jarvis.analysis.thermodynamics.energetics import (
    form_enp,
    get_twod_defect_energy,
    PhaseDiagram,
)
from jarvis.db.figshare import data
import os
from jarvis.io.vasp.outputs import Vasprun

tmp_xml = os.path.join(os.path.dirname(__file__), "JVASP-667_C_C_c.xml")
vrun = Vasprun(tmp_xml)


def test_get_twod_defect_energy():
    Ef = get_twod_defect_energy(vrun=vrun, jid="JVASP-667", atom="C")
    print(Ef)


def test_form_enp():
    atoms = vrun.all_structures[-1]
    total_energy = vrun.final_energy
    Ef = form_enp(atoms=atoms, total_energy=total_energy)
    print(Ef)


def phasediag():
    system = ["Al", "O"]
    system = ["Bi", "Se"]
    system = ["Cu", "Au"]
    system = ["Al", "Ni"]
    system = ["Ga", "Al", "N"]
    system = ["Ni", "Al", "O"]
    system = ["Ni", "Fe", "Cr"]
    system = ["O", "Al", "Ni", "Cu"]
    system=['Mo','Se']
    dft_3d=data('dft_3d')
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
        x, verbose=False, only_plot_stable=True, only_label_stable=True
    )
    for i in x:
        energy, indices, coefs = pd.decompose(i[0])

    pd.plot()


# test_get_twod_defect_energy()
# test_form_enp()
