from jarvis.db.figshare import data, get_jid_data

from jarvis.core.atoms import Atoms, crop_square
from jarvis.io.prismatic.inputs import write_prismatic_xyz

a = Atoms.from_dict(get_jid_data("JVASP-667")["atoms"])


def test_inputs():
    c = crop_square(a)
    lines = write_prismatic_xyz(c)
