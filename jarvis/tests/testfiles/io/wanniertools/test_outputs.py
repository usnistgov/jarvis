import os
from jarvis.io.wanniertools.outputs import WTOut, parse_nodes_dat, parse_chern_dat

wt_out = os.path.join(os.path.dirname(__file__), "WT.out")
wt_nodes = os.path.join(os.path.dirname(__file__), "Nodes.dat")
chrn1 = os.path.join(os.path.dirname(__file__), "wanniercenter3D_Chern.dat")


def test_output():
    z2 = WTOut(path=wt_out).get_z2_index()
    chrn = WTOut(path=wt_out).get_chern_number()
    parse_nodes_dat(wt_nodes)
    assert (z2, chrn[0]) == ("1;0,0,0", 0.0)
    x = parse_chern_dat(chern_dat=chrn1)

#    print(z2)
#    print(chrn)
# test_output()
