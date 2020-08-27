import numpy as np
def qpoint(force_constant=[],qpt=[.1,.1,.1]):
    # FC as natons x natons x 3 x3
    qpt=np.array(qpt)
    exp_iqpt = np.exp(1.0j*qpt)
    dmat=force_constant*exp_iqpt
    val, vect = np.linalg.eigh(dmat)
    return vals,vects
