from jarvis.analysis.stm.tersoff_hamann import TersoffHamannSTM
import matplotlib.pyplot as plt
import os
name = os.path.join(os.path.dirname(__file__), "PARCHG")


def test_th_stm():
    plt.switch_backend("agg")
    TH_STM = TersoffHamannSTM(chg_name=name) 
    t1 = TH_STM.constant_height()
    t2 = TH_STM.constant_current()
