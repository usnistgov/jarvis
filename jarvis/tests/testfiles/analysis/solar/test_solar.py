from jarvis.analysis.solarefficiency.solar import SolarEfficiency
from jarvis.io.vasp.outputs import Vasprun
import os
import tarfile
import matplotlib.pyplot as plt
plt.switch_backend('agg')

example_fold_tgz = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88.tgz",
)


example_fold = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
)

if not os.path.isdir(example_fold):
    tar = tarfile.open(example_fold_tgz)
    tar.extractall(example_fold)
    tar.close()




vrun_path = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
    "SiOptb88",
    "MAIN-MBJ-bulk@mp_149",
    "vasprun.xml",
)


def test_solar():
    v = Vasprun(vrun_path)
    dirgap = v.get_dir_gap
    indirgap = v.get_indir_gap[0]
    en, abz = v.avg_absorption_coefficient
    abz = abz * 100
    eff_slme = SolarEfficiency().slme(
        en, abz, indirgap, indirgap, plot_current_voltage=True
    )
    # print("SLME", 100 * eff)
    eff_sq = SolarEfficiency().calculate_SQ(indirgap,  plot_current_voltage=True)
    # print("SQ", 100 * eff)
    # assert (round(100 * eff_slme, 2), round(100 * eff_sq, 2)) == (33.23, 32.93)
    cmd='rm sq.png slme.png'
    os.system(cmd)

# test_solar()
