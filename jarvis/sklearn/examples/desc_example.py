# coding: utf-8
from jarvis.sklearn.get_desc import get_comp_descp
from pymatgen.core.structure import Structure
from monty.json import MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline
import numpy as np
import sys, os
from jarvis.db.static.explore_db import get_3d_dataset, get_2d_dataset

pos = os.path.join(
    os.path.dirname(__file__), "..", "..", "vasp", "examples", "SiOptb88", "POSCAR"
)


# STEP-1: Check if the descriptor array length for POSCAR is 1557 as
# mentioned in https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.083801
s = Structure.from_file(pos)
X = get_comp_descp(s)
print(X, len(X))
# 1557

# STEP-2: Download JARVIS-DFT database for bulk materials
# from https://figshare.com/articles/jdft_3d-7-7-2018_json/6815699

# STEP-3: load data
data = get_3d_dataset()  # loadfn("jdft_3d-7-7-2018.json", cls=MontyDecoder)
print("data length", len(data))
print("available keys", data[0].keys())


# Simple function to check if the data is available/float
def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


# STEP-4:  Let's train a Formation energy ML
# for i in data:
#  print (i['form_enp'],i['mpid'],i['jid'],i['final_str'],type(i['final_str']))
#  sys.exit()

# Making descriptors and target data exportable to sklearn
X = []
Y = []
for ii, i in enumerate(data):
    if ii < 5:
        y = i["form_enp"]
        if isfloat(y):
            x = get_comp_descp(i["final_str"])
            X.append(x)
            Y.append(y)

# STEP-5: Fit and predict
X = np.array(X).astype(np.float)
Y = np.array(Y).astype(np.float)
est = GradientBoostingRegressor()
pipe = Pipeline([("fs", VarianceThreshold()), ("est", est)])
pipe.fit(X, Y)

test_x = get_comp_descp(data[-1]["final_str"])
print(pipe)
print(pipe.predict([test_x]), data[-1]["form_enp"])
