import numpy as np
from jarvis.db.webpages import Webpage

# Examples of JARVIS IDs with computed Raman Activities
jids_raman_dat = [
    "JVASP-12038",
    "JVASP-75",
    "JVASP-1002",
]


def clean_up_data(dat={}):
    x = {}
    for i, j in dat.items():
        x[i] = np.array(j.strip('"').split(","), dtype="float").tolist()
    return x


mem = []
x = []
# Activities to Intensities
# *************************
# J. At. Mol. Sci. 3 (2012) 1-22
# Intensity = const0 * (nu0-nu)**4/nu * 1/(1-exp(-h*c*nu/(kB*T))) * activity
# Data for computing intensities, as in Materials Project:
#      T=300K,
#      nu0=18796.9925 cm-1
temp = 300.0
nu0 = 18796.9925
# hc = Plank's constant * speed of light
# hc_over_kB = h*c/Boltzman_constant
hc_over_kB = 1.4388003592  # cm/s
const1 = hc_over_kB / temp
# A good, aritrary choice for const0 is 10**(-14)
const0 = 10 ** (-14)
for ii, i in enumerate(jids_raman_dat):
    w = Webpage(jid=i)
    dat = w.data["basic_info"]["raman_dat"]
    if dat is not None:
        info = {}
        info["jid"] = i
        info["dat"] = clean_up_data(dat)
        mem.append(info)
        n_modes = len(info["dat"]["frequencies"])
        print(ii, i, len(mem), n_modes)
        print("wave#[cm-1]    activity    intensity")
        for kk in range(n_modes):
            #print(info["dat"]["frequencies"][kk], info["dat"]["indices"][kk])
            # Activities to Intensities
            activity = info["dat"]["indices"][kk]
            nu = info["dat"]["frequencies"][kk]
            AA1 = const0 * (nu0 - nu) ** 4 / nu
            exp1 = np.exp(-hc_over_kB * nu / temp)
            AA2 = 1.0 / (1.0 - exp1)
            intensity = AA1 * AA2 * activity
            print(nu, activity, intensity)
        print("")
