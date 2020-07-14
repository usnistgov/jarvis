from jarvis.analysis.stm.tersoff_hamann import TersoffHamannSTM
import matplotlib.pyplot as plt
import os
name = os.path.join(os.path.dirname(__file__), "PARCHG")


def test_th_stm():
    plt.switch_backend("agg")
    TH_STM1 = TersoffHamannSTM(chg_name=name) 
    t1 = TH_STM1.constant_height()
    TH_STM2 = TersoffHamannSTM(chg_name=name) 
    t2 = TH_STM2.constant_current()
