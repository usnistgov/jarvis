from jarvis.db.figshare import data, get_jid_data
import matplotlib.pyplot as plt
from jarvis.core.atoms import Atoms, crop_square

plt.switch_backend("agg")
from jarvis.analysis.stem.convolution_apprx import STEMConv


def test_conv_method():

    plt.switch_backend("agg")

    a = Atoms.from_dict(get_jid_data("JVASP-667")["atoms"])
    c = crop_square(a)
    # c = a.make_supercell_matrix([2, 2, 1])
    p = STEMConv(atoms=c).simulate_surface(c)

    # plt.imshow(p, interpolation="gaussian", cmap="plasma")
    # plt.savefig("stem_example.png")
    # plt.close()


# test_conv_method()
