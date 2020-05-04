import os
from jarvis.io.wanniertools.outputs import WTOut

wt_out = os.path.join(os.path.dirname(__file__), "WT.out")


def test_output():
    z2 = WTOut(path=wt_out).get_z2_index()
    chrn = WTOut(path=wt_out).get_chern_number()
    assert (z2, chrn[0]) == ("1;0,0,0", 0.0)


#    print(z2)
#    print(chrn)
# test_output()
