from jarvis.analysis.periodic.ptable import plot_ptable_trend
import os


def test_make_bokeh():
    plot_ptable_trend(save_plot=True)
    plot_ptable_trend(save_plot=False)
    os.remove("ptable.html")
