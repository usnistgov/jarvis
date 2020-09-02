"""Module for post-processing phonopy outputs."""
import yaml
import matplotlib.pyplot as plt
from yaml import Loader
import numpy as np


def bandstructure_plot(band_yaml="", plot=False):
    """Get phonopy bandstructure info."""
    with open(band_yaml, "r") as f:
        data = yaml.load(f, Loader=Loader)
    frequencies = []
    distances = []
    qpoints = []
    label_points = []
    labels = []
    for j, v in enumerate(data["phonon"]):
        if "label" in v and v["label"] != "None":
            labels.append(v["label"])
            label_points.append(v["distance"])
        frequencies.append([f["frequency"] for f in v["band"]])
        qpoints.append(v["q-position"])
        distances.append(v["distance"])
    if plot:
        for i in range(np.array(frequencies).shape[1]):
            plt.plot(distances, np.array(frequencies)[:, i])
        plt.xticks(label_points, labels)
    return frequencies, distances, labels, label_points


def total_dos(tot_dos="", plot=False):
    """Get total dos info."""
    f = open(tot_dos, "r")
    freq = []
    pdos = []
    for lines in f.readlines():
        if not str(lines.split()[0]).startswith("#"):
            #   print (lines)
            # else:
            freq.append(float(lines.split()[0]))
            pdos.append(float(lines.split()[1]))
    if plot:
        plt.plot(freq, pdos)
    return freq, pdos
