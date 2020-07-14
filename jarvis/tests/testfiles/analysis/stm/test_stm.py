from jarvis.analysis.stm.tersoff_hamann import TersoffHamannSTM
import matplotlib.pyplot as plt
import os
f = os.path.join(os.path.dirname(__file__), "PARCHG")


def test_th_stm():
    plt.switch_backend("agg")
    t = TersoffHamannSTM(chg_name=f).constant_height()
    t = TersoffHamannSTM(chg_name=f).constant_current()
