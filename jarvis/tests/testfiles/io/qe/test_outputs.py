import os
from jarvis.io.qe.outputs import QEout, DataFileSchema, ProjHamXml
import matplotlib.pyplot as plt

plt.switch_backend("agg")
out = os.path.join(os.path.dirname(__file__), "qe.out")
xml = os.path.join(os.path.dirname(__file__), "data-file-schema.xml")
projham_xml = os.path.join(os.path.dirname(__file__), "projham_K.xml")


def test_outputs():
    en = QEout(filename=out).get_total_energy()
    print((en))
    td = QEout(filename=out).to_dict()
    fd = QEout.from_dict(td)
    # assert en == -19.11812163
    en = QEout(filename=out).get_band_enegies()
    print((en), len(en))
    # assert en[0][0] == -5.8325
    en = QEout(filename=out).get_efermi()
    print((en))
    # assert en == 6.4236
    dxml = DataFileSchema(filename=xml)
    print("final energy", dxml.final_energy)
    print("final energy breakdown", dxml.final_energy_breakdown)
    print("forces=", dxml.forces)
    print("stress", dxml.stress)
    print("magnetization", dxml.magnetization)
    print(dxml.indir_gap)
    print(dxml.functional)
    print(dxml.initial_structure)
    print(dxml.initial_structure)
    print(dxml.final_structure)
    print(dxml.bandstruct_eigvals(plot=True))
    cmd = "rm band.png"
    projham = ProjHamXml(filename=projham_xml).get_tight_binding()
    energies, dos, pdos, names = ProjHamXml(filename=projham_xml).dos()
    print("energies", energies)
    os.system(cmd)


# test_outputs()
