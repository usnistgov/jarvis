"""Helper functions for ML applications."""

from jarvis.db.figshare import data
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
import numpy as np

typical_data_ranges = {
    "formation_energy_peratom": [-5, 5],
    "optb88vdw_bandgap": [0, 10],
    "mbj_bandgap": [0, 10],
    "bulk_modulus_kv": [0, 250],
    "shear_modulus_gv": [0, 250],
    "epsx": [0, 60],
    "epsy": [0, 60],
    "epsz": [0, 60],
    "mepsx": [0, 60],
    "mepsy": [0, 60],
    "mepsz": [0, 60],
    "n-Seebeck": [-600, 10],
    "n-powerfact": [0, 5000],
    "p-Seebeck": [-10, 600],
    "p-powerfact": [0, 5000],
    "slme": [0, 40],
    "spillage": [0, 4],
    "encut": [0, 2000],
    "kpoint_length_unit": [0, 200],
    "dfpt_piezo_max_dielectric": [0, 100],
    "dfpt_piezo_max_dij": [0, 3000],
    "dfpt_piezo_max_eij": [0, 10],
    "ehull": [0, 1],
    "electron_avg_effective_masses_300K": [0, 3],
    "hole_avg_effective_masses_300K": [0, 3],
    "exfoliation_energy": [0, 1000],
    "magmom_oszicar": [0, 10],
    "max_ir_mode": [0, 4000],
    "total_energy_per_atom": [-10, 3],
}


def get_ml_data(
    ml_property="formation_energy_peratom",
    dataset="cfid_3d",
    data_ranges=typical_data_ranges,
):
    """
    Provide arrays/pandas-dataframe as input for ML algorithms.

    Args:

        ml_property: target property to train

        data_ranges: range for filtering data

        dataset: dataset available in jarvis or other array

    Returns:
           X, Y , ids
    """
    import pandas as pd

    if isinstance(dataset, str):
        dataml = data(dataset)
        df = pd.DataFrame(dataml)
    else:
        df = pd.DataFrame(dataset)

    x = []
    y = []
    jid = []
    df2 = df[["desc", "jid", ml_property]].replace("na", np.nan).dropna()
    for ii, i in df2.iterrows():
        if data_ranges is None:
            if (
                len(i["desc"]) == 1557
                and float(i[ml_property]) != float("inf")
                and i[ml_property] != "na"
            ):
                x.append(i["desc"])
                y.append(i[ml_property])
                jid.append(i["jid"])
        else:
            if (
                len(i["desc"]) == 1557
                and float(i[ml_property]) != float("inf")
                and i[ml_property] != "na"
                and float(i[ml_property]) <= data_ranges[ml_property][1]
                and float(i[ml_property]) >= data_ranges[ml_property][0]
            ):
                x.append(i["desc"])
                y.append(i[ml_property])
                jid.append(i["jid"])
    return np.array(x, dtype="float"), np.array(y, dtype="float"), jid


def mean_absolute_deviation(data, axis=None):
    """Get Mean absolute deviation."""
    return np.mean(np.absolute(data - np.mean(np.array(data), axis)), axis)


def regr_scores(test, pred):
    """Provide generic regresion scores.

    Args:

       pred: predicted values

       test: held data for testing

    Returns:
        info: with metrics
    """
    rmse = np.sqrt(mean_squared_error(test, pred))
    r2 = r2_score(test, pred)
    mae = mean_absolute_error(test, pred)
    info = {}
    info["mae"] = mae
    info["rmse"] = rmse
    info["r2"] = r2
    info["test"] = test
    info["pred"] = pred
    return info


def binary_class_dat(X=[], Y=[], tol=0.1):
    """
    Categorize a continous dataset in 1/0 with a threshold "tol".

    TODO: replace with OneHotEncoder
    """
    Y1 = []
    for i, j in zip(X, Y):
        if j >= tol:
            Y1.append(1)
        else:
            Y1.append(0)
    return X, Y1
