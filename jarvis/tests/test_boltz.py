import os
from jarvis.boltztrap.boltztrap import get_prop,boltz_run
mainrun = os.path.join(
    os.path.dirname(__file__),
    "..",
    "vasp",
    "examples",
    "SiOptb88",
    "MAIN-RELAX-bulk@mp_149",
    "vasprun.xml",
)
"""
def test_boltz1():
    b = boltz_run(mainrun)
    val = get_prop(b, prop="zt")
    assert val[0]==1.5913528701329165

"""
