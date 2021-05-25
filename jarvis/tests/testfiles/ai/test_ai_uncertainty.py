from jarvis.ai.uncertainty import lgbm_quantile_uncertainty as uq
from jarvis.ai.pkgs.utils import get_ml_data
import os

# STEP-1: Getting Data
# ********************
property = "exfoliation_energy"

x, y, jid = get_ml_data(dataset="cfid_3d", ml_property=property)

x = x[0:100]
y = y[0:100]
jid = jid[0:100]

# STEP-2: Quantile regression to make predictions and determine
#         prediction intervals for predicted data
# **************************************************************

# Search parameters
scoring = "neg_mean_absolute_error"
cv = 2
n_jobs = -1
# n_iter = how many parameter combination to try in the search
n_iter = 10
random_state = 508842607

info = uq.quantile_regr_predint(
    x, y, jid, cv, n_jobs, n_iter, random_state, scoring, property
)

# TEST
assert info["MAE_Lower"] < 200.0
assert info["MAE_Mid"] < 200.0
assert info["MAE_Upper"] < 200.0
assert info["MAE_Error"] < 200.0
cmd ='rm *.pk *.dat'
os.system(cmd)
