import os
from jarvis.io.qe.outputs import QEout, DataFileSchema
import matplotlib.pyplot as plt

plt.switch_backend("agg")
out = os.path.join(os.path.dirname(__file__), "qe.out")
xml = os.path.join(os.path.dirname(__file__), "data-file-schema.xml")


def test_outputs():
    en = QEout(filename=out).get_total_energy()
    print((en))
    td = QEout(filename=out).to_dict()
    fd = QEout.from_dict(td)
    assert en == -19.11812163
    en = QEout(filename=out).get_band_enegies()
    print((en), len(en))
    assert en[0][0] == -5.8325
    en = QEout(filename=out).get_efermi()
    print((en))
    assert en == 6.4236
    dxml = DataFileSchema(filename=xml)
    print(dxml.final_energy)
    print(dxml.forces)
    print(dxml.final_structure)
    print(dxml.bandstruct_eigvals(plot=True))
    cmd = "rm band.png"
    os.system(cmd)
    