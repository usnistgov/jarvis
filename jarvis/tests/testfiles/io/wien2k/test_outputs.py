from jarvis.io.wien2k.outputs import band_eigvals, read_spaghetti_ene, read_scf, read_band_energy
import os
import matplotlib.pyplot as plt
plt.switch_backend('agg')
energy_file = os.path.join(os.path.dirname(__file__), "FeSe.energy")
scf_file = os.path.join(os.path.dirname(__file__), "FeSe.scf")
ene_file = os.path.join(os.path.dirname(__file__), "ICSD-42662.spaghetti_ene")


def test_out():
    eigs = read_band_energy(energy_file=energy_file)
    info = read_scf(scf_file=scf_file)
    band_eigvals(energy_file=energy_file, plot=True)
    a,b=read_spaghetti_ene(ene_file)
