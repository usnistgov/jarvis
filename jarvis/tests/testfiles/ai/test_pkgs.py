from jarvis.ai.pkgs.sklearn.regression import regression
from jarvis.ai.pkgs.sklearn.classification import classification
from jarvis.ai.pkgs.sklearn.hyper_params import (
    classification_regression_params,
)
from jarvis.ai.pkgs.utils import get_ml_data, binary_class_dat
from jarvis.ai.pkgs.lgbm.regression import regression as l_regression
from jarvis.ai.pkgs.lgbm.regression import parameters_dict as l_params
from jarvis.ai.pkgs.lgbm.classification import (
    classification as l_classification,
)
from jarvis.ai.descriptors.cfid import feat_names
from lightgbm import LGBMClassifier
import matplotlib.pyplot as plt
import os

plt.switch_backend("agg")

property = "exfoliation_energy"

X, Y, jid = get_ml_data(dataset="cfid_3d", ml_property=property)

X = X[0:100]
Y = Y[0:100]
jid = jid[0:100]


def test_lgbm_regression():
    params = l_params()[property]
    names = feat_names()
    info = l_regression(
        X=X,
        Y=Y,
        jid=jid,
        config=params,
        feat_names=names,
        plot=True,
        save_model=True,
    )
    assert info["reg_scores"]["mae"] < 200.0
    info = l_regression(
        X=X,
        Y=Y,
        jid=jid,
        config=params,
        preprocess=False,
        feat_names=names,
        plot=True,
        save_model=True,
    )
    cmd = 'rm my_model* *.png'
    os.system(cmd)

def test_lgbm_classification():
    property = "exfoliation_energy"
    tol = 100
    models = [LGBMClassifier(n_estimators=10, num_leaves=2)]
    info = l_classification(
        X=X,
        Y=Y,
        plot=True,
        models=models,
        preprocess=True,
        save_model=True,
        tol=tol,
    )
    cmd = 'rm my_model* *.png'
    os.system(cmd)
    assert info["LGBMClassifier"]["roc_auc"][0] > 0.0


def test_hyperparams():
    x, y = classification_regression_params()


def test_sklearn_simple_regression():

    info = regression(X, Y, plot=True)
    assert info["GradientBoostingRegressor"]["mae"] < 200.0


def test_sklearn_simple_classification():
    X_class, Y_class = binary_class_dat(X=X, Y=Y, tol=100)
    info = classification(X=X, Y=Y, tol=100, plot=True)
    assert (info["GradientBoostingClassifier"]["roc_auc"][0]) > 0.0
