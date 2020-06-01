from jarvis.ai.pkgs.sklearn.regression import regression
from jarvis.ai.pkgs.sklearn.classification import  classification
from jarvis.ai.pkgs.sklearn.hyper_params import classification_regression_params
from jarvis.ai.pkgs.utils import get_ml_data, binary_class_dat
from jarvis.ai.pkgs.lgbm.regression import regression as l_regression
from jarvis.ai.pkgs.lgbm.regression import parameters_dict as l_params
from jarvis.ai.pkgs.lgbm.classification import classification as l_classification
from jarvis.ai.descriptors.cfid import feat_names
from lightgbm import LGBMClassifier

property = "exfoliation_energy"

X, Y, jid  = get_ml_data(dataset = 'cfid_3d', ml_property=property)


def test_lgbm_regression():
    params = l_params()[property]
    names = feat_names()
    info = l_regression(X=X, Y=Y,jid=jid, config=params, feat_names=names)
    assert info['reg_scores']['mae']<50.0

"""    
def test_lgbm_classification():
    property = "exfoliation_energy"
    tol=100
    models = [LGBMClassifier()]
    info = l_classification(X=X, Y=Y, models=models, preprocess=True, save_model=False, tol=tol)
    assert info['LGBMClassifier']['roc_auc'][0]>0.5
"""

def test_hyperparams():
   x,y  = classification_regression_params()


def test_sklearn_simple_regression():

    info = regression(X, Y)
    assert info['GradientBoostingRegressor']['mae']<50.0

def test_sklearn_simple_classification():
    X_class, Y_class = binary_class_dat(X=X, Y=Y, tol=100)
    info = classification(X=X,Y=Y,tol=100)
    assert (info['GradientBoostingClassifier']['roc_auc'][0])>0.75



