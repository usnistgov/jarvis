from jarvis.analysis.periodic.ptable import plot_ptable_trend
import os

filename = (os.path.join(os.path.dirname(__file__), "born.csv"))
def test_make_bokeh():
    plot_ptable_trend(save_plot=True)
    plot_ptable_trend(save_plot=True,log_scale=1)
    plot_ptable_trend(save_plot=False, input_file=filename)
    os.remove("ptable.html")
    