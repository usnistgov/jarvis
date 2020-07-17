import os
from jarvis.io.qe.outputs import QEout

out = os.path.join(os.path.dirname(__file__), "qe.out")


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
