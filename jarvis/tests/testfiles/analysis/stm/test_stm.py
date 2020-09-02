from jarvis.analysis.stm.tersoff_hamann import TersoffHamannSTM
import matplotlib.pyplot as plt
import os
name = os.path.join(os.path.dirname(__file__), "PARCHG")

from io import BytesIO
def test_th_stm():
    plt.switch_backend("agg")
    TH_STM1 = TersoffHamannSTM(chg_name=name)
    byte_io = BytesIO() 
    t1 = TH_STM1.constant_height(filename=byte_io)
    p=byte_io.getvalue()#.decode('UTF-8')
    #print ('p',p)
    t1 = TH_STM1.constant_height()
    TH_STM2 = TersoffHamannSTM(chg_name=name) 
    t2 = TH_STM2.constant_current()
    t2 = TH_STM2.constant_current(pc=5)
    cmd = 'rm *.png'
    os.system(cmd)

#test_th_stm()
