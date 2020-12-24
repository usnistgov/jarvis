import os
from jarvis.io.wanniertools.outputs import WTOut, parse_nodes_dat

wt_out = os.path.join(os.path.dirname(__file__), "WT.out")
wt_nodes = os.path.join(os.path.dirname(__file__), "Nodes.dat")


def test_output():
    z2 = WTOut(path=wt_out).get_z2_index()
    chrn = WTOut(path=wt_out).get_chern_number()
    parse_nodes_dat(wt_nodes)
    assert (z2, chrn[0]) == ("1;0,0,0", 0.0)


#    print(z2)
#    print(chrn)
# test_output()
