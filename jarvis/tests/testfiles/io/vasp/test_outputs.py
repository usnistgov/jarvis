from jarvis.io.vasp.outputs import Vasprun, Oszicar, Wavecar, Waveder, Chgcar, Outcar
import numpy as np
import os

vrun = Vasprun(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "MAIN-RELAX-bulk@mp_149",
        "vasprun.xml",
    )
)
chg = Chgcar(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "MAIN-RELAX-bulk@mp_149",
        "CHGCAR",
    )
)
loc = Chgcar(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "MAIN-RELAX-bulk@mp_149",
        "LOCPOT",
    )
)
osz = Oszicar(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "MAIN-RELAX-bulk@mp_149",
        "OSZICAR",
    )
)
out = Outcar(
    filename=os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "MAIN-ELASTIC-bulk@mp_149",
        "OUTCAR",
    )
)
wder = Waveder(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "..",
        "examples",
        "vasp",
        "SiOptb88",
        "MAIN-OPTICS-bulk@mp_149",
        "WAVEDER",
    )
)
wf_noso = Wavecar(
    filename=os.path.join(
        os.path.dirname(__file__), "..", "..", "analysis", "spillage", "WAVECAR.nosoc"
    ),
    lsorbit=False,
)
wf_so = Wavecar(
    filename=os.path.join(
        os.path.dirname(__file__), "..", "..", "analysis", "spillage", "WAVECAR.soc"
    ),
    lsorbit=True,
)


def test_chgcar():
    # print (chg.is_spin_polarized(), chg.is_spin_orbit(), np.array(chg.chg).shape)
    assert (chg.is_spin_polarized(), chg.is_spin_orbit(), np.array(chg.chg).shape) == (
        True,
        False,
        (4, 56, 56, 56),
    )


def test_vrun():
    # print ('gapp',round(vrun.get_indir_gap,2))
    assert (round(vrun.get_indir_gap, 2)) == (0.73)


def test_osz():
    assert (float(osz.magnetic_moment)) == (0.0)


def test_out():
    assert (round(out.elastic_props()["KV"], 2)) == (87.27)
