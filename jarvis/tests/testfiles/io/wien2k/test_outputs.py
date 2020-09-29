from jarvis.io.wien2k.outputs import read_scf, read_band_energy
import os

energy_file = os.path.join(os.path.dirname(__file__), "FeSe.energy")
scf_file = os.path.join(os.path.dirname(__file__), "FeSe.scf")


def test_out():
    eigs = read_band_energy(energy_file=energy_file)
    info = read_scf(scf_file=scf_file)
