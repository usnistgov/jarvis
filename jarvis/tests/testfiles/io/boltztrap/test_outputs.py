from jarvis.io.boltztrap.outputs import BoltzTrapOutput
import os

bpath = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
    "MAIN-RELAX-bulk@mp_149",
    "boltztrap",
)


def test_out():
    b = BoltzTrapOutput(path=bpath)
    x = b.to_dict()["condtens_fixdoping"]
    tmp = x["p"][300.0][4.1e-07]["cond"][0]
    assert tmp == 609055450000000.0


# test_out()
