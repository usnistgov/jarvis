from jarvis.core.atoms import crop_squre
from jarvis.db.figshare import data, get_jid_data
import matplotlib.pyplot as plt
from jarvis.core.atoms import Atoms, ase_to_atoms, get_supercell_dims

plt.switch_backend("agg")
from jarvis.analysis.stem.convolution_apprx import (
    simulate_surface,
    probe_profile,
)


def test_conv_method():
    a = Atoms.from_dict(get_jid_data("JVASP-667")["atoms"])
    c = crop_squre(a)
    p = simulate_surface(
        c,
        probe_profile=probe_profile(
            gaussian_width=1, lorentzian_width=1, intensity_ratio=1, n=100
        ),
        shape=[50, 50],
    )
    # plt.imshow(p, interpolation="gaussian", cmap="plasma")
    # plt.savefig("stem_example.png")
    # plt.close()


# test_conv_method()
