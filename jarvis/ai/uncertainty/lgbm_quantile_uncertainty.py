# form_enp,op_gap,mbj_gap,kp_leng,encut,fin_enp,magmom,epsx,epsy,epsz,mepsx,mepsy,mepsz,exfoliation_en,type,mpid,jid

# from jarvis.sklearn.get_desc import get_comp_descp
from monty.serialization import loadfn, MontyEncoder, MontyDecoder, dumpfn
import numpy as np
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
import lightgbm as lgb

# import matplotlib.pyplot as plt
# plt.switch_backend('agg')
import pandas as pd
from sklearn.datasets import load_boston
from sklearn.model_selection import (
    train_test_split,
    learning_curve,
    cross_val_score,
    cross_val_predict,
    GridSearchCV,
    RandomizedSearchCV,
)
import scipy as sp
import time, os, json, pprint
from sklearn.feature_selection import (
    SelectKBest,
    f_classif,
    SelectFromModel,
    VarianceThreshold,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import pickle

# import joblib
import numpy as np
from sklearn.cluster import KMeans
import pandas as pd

# Modeling
from sklearn.base import BaseEstimator
from sklearn.ensemble import GradientBoostingRegressor

# File finding
import glob


# ---------------------------------------------------------------------


def get_number_formula_unit(s=""):

    orig_formula = s.composition.as_dict()
    prim_formula = s.get_primitive_structure().composition.as_dict()
    num_unit = float(orig_formula.values()[0]) / float(
        prim_formula.values()[0]
    )
    return num_unit


def isfloat(value):
    try:
        float(value)
        return True
    except:
        return False
        pass


def jdata(prop=""):
    # d3=loadfn('/rk2/ftavazza/ML/NEW/Quantile/New_quantities/jml_3d-4-26-2020.json',cls=MontyDecoder)
    d3 = loadfn(
        "/users/ftavazza/ML/ruth_New/jml_3d-4-26-2020.json", cls=MontyDecoder
    )
    # limits={'form_enp':[-5,5],'exfoliation_en':[0,1000],'el_mass_x':[0,100],'el_mass_y':[0,100],'el_mass_z':[0,100],'hl_mass_x':[0,100],'hl_mass_y':[0,100],'hl_mass_z':[0,100],'magmom':[0,10],'fin_en':[0,0],'op_gap':[-1,10],'kv':[-1,500],'gv':[-15,510],'encut':[500,3000],'kp_leng':[40,200],'epsx':[0,50],'epsy':[0,50],'epsz':[0,50],'mepsx':[0,50],'mepsy':[0,50],'mepsz':[0,50],'mbj_gap':[-1,10]}

    limits = {
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
    X = []
    Y = []
    jid = []
    for ii, i in enumerate(d3):
        # if ii<100:
        y = i[prop]
        if isfloat(y):
            y = float(y)
            if y >= limits[prop][0] and y <= limits[prop][1]:
                x = i["desc"]  # get_comp_descp(i['final_str'])
                # if len(x)==1557 and any(np.isnan(x) for x in x.flatten())==False:  # and np.isfinite(y):
                if (
                    len(x) == 1557
                ):  # and any(np.isnan(x) for x in x.flatten())==False:  # and np.isfinite(y):
                    if "eps" in prop:
                        y = np.sqrt(float(y))
                    if "mag" in prop:
                        num = get_number_formula_unit(i["strt"])
                        y = float(abs(y)) / float(num)

                    X.append(x)
                    Y.append(y)
                    jid.append(i["jid"])
    print("max,min", max(Y), min(Y))
    print("Prop=", prop, len(X), len(Y))
    X = np.array(X).astype(np.float64)
    Y = np.array(Y).astype(np.float64)
    return X, Y, jid


from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.ensemble import (
    RandomForestRegressor,
    GradientBoostingRegressor,
    AdaBoostRegressor,
)
from sklearn import datasets, svm
from sklearn.kernel_approximation import Nystroem
from sklearn.svm import SVR, LinearSVR, SVC, LinearSVC
from sklearn.linear_model import Lasso, LinearRegression, LogisticRegression
from sklearn.linear_model import LassoCV, LassoLarsCV, LogisticRegressionCV
from sklearn.linear_model import RidgeCV
from sklearn.kernel_ridge import KernelRidge
from sklearn.neural_network import MLPRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.multiclass import OneVsRestClassifier
from sklearn.datasets import load_boston
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline
from sklearn.neural_network import MLPClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.ensemble import (
    RandomForestClassifier,
    AdaBoostClassifier,
    GradientBoostingClassifier,
)
from sklearn.svm import SVC


# http://scikit-learn.org/stable/modules/classes.html#module-sklearn.model_selection


def cluster_cv(X=[], y=[], ids=[], technique=KMeans()):
    labels = technique.fit_predict(X)
    clust_labels = np.unique(labels)
    last_clust = np.unique(labels)[-1]
    test_x = []
    test_y = []
    test_z = []
    train_x = []
    train_y = []
    train_z = []

    for i, j, k, l in zip(labels, X, y, ids):
        if i == last_clust:
            test_x.append(j)
            test_y.append(k)
            test_z.append(l)
        else:
            train_x.append(j)
            train_y.append(k)
            train_z.append(l)

    return (
        np.array(train_x),
        np.array(train_y),
        np.array(train_z),
        np.array(test_x),
        np.array(test_y),
        np.array(test_z),
        labels,
    )


def plot_learning_curve(
    estimator=None,
    scoring="neg_mean_absolute_error",
    title="",
    X=None,
    y=None,
    ylim=None,
    cv=5,
    # n_jobs=-1, train_sizes=np.linspace(.1, 1.0, 10),fname='fig.png'):
    n_jobs=-1,
    train_sizes=np.linspace(0.01, 1.0, 50),
    fname="fig.png",
):
    """
    Taken from scikit-learn, added ways to store results in JSON format

    """

    plt.figure()
    # fname='fig.png'
    plt.title(title)
    mem = {}
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    train_sizes, train_scores, test_scores = learning_curve(
        estimator,
        X,
        y,
        cv=cv,
        n_jobs=n_jobs,
        train_sizes=train_sizes,
        scoring=scoring,
    )
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    plt.grid()
    mem["train_sizes"] = train_sizes
    mem["train_scores_mean"] = train_scores_mean
    mem["test_scores_mean"] = test_scores_mean
    mem["train_scores_std"] = train_scores_std
    mem["test_scores_std"] = test_scores_std

    plt.fill_between(
        train_sizes,
        train_scores_mean - train_scores_std,
        train_scores_mean + train_scores_std,
        alpha=0.1,
        color="r",
    )
    plt.fill_between(
        train_sizes,
        test_scores_mean - test_scores_std,
        test_scores_mean + test_scores_std,
        alpha=0.1,
        color="g",
    )
    plt.plot(
        train_sizes, train_scores_mean, "o-", color="r", label="Training score"
    )
    plt.plot(
        train_sizes,
        test_scores_mean,
        "o-",
        color="g",
        label="Cross-validation score",
    )

    plt.legend(loc="best")
    fname = str(fname) + str("_learn_1.png")
    plt.show()
    plt.savefig(fname)
    plt.close()

    data = str(fname) + str("learn_data")
    f = open(data, "w")
    line = (
        str("train_sizes ")
        + str("train_scores_mean ")
        + str("train_scores_std ")
        + str("test_scores_mean ")
        + str("test_scores_std ")
        + "\n"
    )
    f.write(line)
    for i, j, k, l, m in zip(
        train_sizes,
        train_scores_mean,
        train_scores_std,
        test_scores_mean,
        test_scores_std,
    ):
        line = (
            str(i)
            + str(" ")
            + str(j)
            + str(" ")
            + str(k)
            + str(" ")
            + str(l)
            + str(" ")
            + str(m)
            + str("\n")
        )
        f.write(line)
    f.close()

    plt.figure()
    plt.title(title)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Training examples")
    plt.ylabel("Score")
    plt.grid()
    plt.fill_between(
        train_sizes,
        test_scores_mean - test_scores_std,
        test_scores_mean + test_scores_std,
        alpha=0.1,
        color="g",
    )
    plt.plot(
        train_sizes,
        test_scores_mean,
        "o-",
        color="g",
        label="Cross-validation score",
    )
    plt.legend(loc="best")
    fname = str(fname) + str("_learn_2.png")
    plt.show()
    plt.savefig(fname)
    plt.close()
    return mem


# def regr_scores(pred,test):
def regr_scores(test, pred):

    """
    Generic regresion scores

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


def get_lgbm(
    train_x,
    val_x,
    train_y,
    val_y,
    cv,
    n_jobs,
    scoring,
    n_iter,
    objective,
    alpha,
    random_state,
):
    # def get_lgbm(train_x, val_x, train_y,val_y,cv,n_jobs,scoring,n_iter):
    """
    Train a lightgbm model

    Args:
        train_x: samples used for trainiing
        val_x: validation set
        train_y: train targets
        val_y: validation targets
        cv: # of cross-validations
        n_jobs: for making the job parallel
        scoring: scoring function to use such as MAE
    Returns:
           Best estimator
    """

    # Get converged boosting iterations with high learning rate, MAE as the convergence crietria
    lgbm = lgb.LGBMRegressor(
        n_estimators=500,
        learning_rate=0.1,
        max_depth=5,
        num_leaves=100,
        objective=objective,
        # min_data_in_leaf=2,
        n_jobs=-1,
        alpha=alpha,
        random_state=random_state,
        verbose=-1,
    )

    lgbm.fit(
        train_x,
        train_y,
        eval_set=[(val_x, val_y)],
        eval_metric="mae",
        # eval_metric='l1',
        early_stopping_rounds=10,
    )
    num_iteration = lgbm.best_iteration_
    print("num_iteration", num_iteration)
    print("in randomsearch cv")
    # Generally thousands of randomized search for optimal parameters
    # learning rate and num_leaves are very important
    param_dist = {
        #'boosting_type': [ 'dart'],
        #'boosting_type': ['gbdt', 'dart', 'rf'],
        #'num_leaves': sp.stats.randint(2, 1001),
        #'subsample_for_bin': sp.stats.randint(10, 1001),
        #'min_split_gain': sp.stats.uniform(0, 5.0),
        #'min_child_weight': sp.stats.uniform(1e-6, 1e-2),
        #'reg_alpha': sp.stats.uniform(0, 1e-2),
        #'reg_lambda': sp.stats.uniform(0, 1e-2),
        #'tree_learner': ['data', 'feature', 'serial', 'voting' ],
        #'application': ['regression_l1', 'regression_l2', 'regression'],
        #'bagging_freq': sp.stats.randint(1, 11),
        #'bagging_fraction': sp.stats.uniform(.1, 0.9),
        #'feature_fraction': sp.stats.uniform(.1, 0.9),
        #'learning_rate': sp.stats.uniform(1e-3, 0.9),
        #'est__num_leaves': [2,8,16],
        #'est__min_data_in_leaf': [1,2,4],
        #'est__learning_rate': [0.005,0.01,0.1],
        #'est__max_depth': [1,3,5], #sp.stats.randint(1, 501),
        #'est__n_estimators': [num_iteration,2*num_iteration,5*num_iteration],#sp.stats.randint(100, 20001),
        #'gpu_use_dp': [True, False],
        "est__min_data_in_leaf": sp.stats.randint(5, 20),
        "est__n_estimators": sp.stats.randint(500, 2000),
        "est__num_leaves": sp.stats.randint(100, 500),
        "est__max_depth": sp.stats.randint(8, 20),
        "est__learning_rate": sp.stats.uniform(5e-3, 0.5),
    }
    lgbm = lgb.LGBMRegressor(
        objective=objective,
        # device='gpu',
        # n_estimators=num_iteration,
        n_jobs=n_jobs,
        alpha=alpha,
        verbose=-1,
    )
    pipe = Pipeline(
        [
            ("stdscal", StandardScaler()),
            ("vart", VarianceThreshold(1e-4)),
            ("est", lgbm),
        ]
    )

    # n_iter=10
    # Increase n_iter
    rscv = RandomizedSearchCV(
        estimator=pipe,
        param_distributions=param_dist,
        cv=cv,
        scoring=scoring,
        n_iter=n_iter,
        n_jobs=n_jobs,
        verbose=3,
        random_state=random_state,
        refit=True,
    )
    model = rscv.fit(train_x, train_y)
    print("Best Score: ", model.best_score_)
    print("Best Params: ", model.best_params_)
    print("Best Estimator: ", model.best_estimator_)
    # return model.best_estimator_
    return model


# Main-function
def run(
    version="version_1",
    scoring="neg_mean_absolute_error",
    cv=5,
    n_jobs=1,
    prop="op_gap",
    do_cv=False,
):
    """
    Generic run function to train-test split,
    find optimum number of boosting operations,
    the hyperparameter optimization,
    learning curves, n-fold cross-validations,
    saving parameters and results


    Args:
        version: user defined name
        scoring: evaluation metric
        cv: # of cross-validation
        n_jobs: for running in parallel
        prop: property to train
        do_cv: whether to perform cross validation

    """
    name = str(version) + str("_") + str(prop)
    # Make a directory for storing results
    dir_name = str(os.getcwd()) + str("/") + str(name)
    # if not os.path.exists(dir_name): os.makedirs(dir_name)
    # os.chdir(dir_name)

    info = {}
    tmp_time = time.time()

    # STEP-1: Data for a particular model
    # ***********************************
    # Toy example with boston data
    # boston=load_boston()
    # x, y = boston['data'], boston['target']

    x, y, jid = jdata(prop)

    # STEP-2: Splitting the data
    # ***************************
    # 90-10% split for train test

    X_train, X_test, y_train, y_test, jid_train, jid_test = train_test_split(
        x, y, jid, random_state=1, test_size=0.1
    )
    # print ('lenx len y',len(x[0]),len(y))

    # STEP-3: Use a specific ML model
    # ********************************
    data_mae = []
    data_rmse = []
    data_r2 = []
    data_test = {}
    data_pred = {}
    model_list = [
        "KernelRidge",
        "MLPRegressor",
        "GaussianProcessRegressor",
        "RandomForestRegressor",
        "GradientBoostingRegressor",
        "AdaBoostRegressor",
        "LinearRegression",
        "DecisionTreeRegressor",
        "Lasso",
        "LinearSVR",
        "LogisticRegression",
    ]
    model_list = [
        "KernelRidge",
        "MLPRegressor",
        "GaussianProcessRegressor",
        "RandomForestRegressor",
        "GradientBoostingRegressor",
        "AdaBoostRegressor",
        "LinearRegression",
        "DecisionTreeRegressor",
        "Lasso",
    ]

    from sklearn import linear_model
    from sklearn import datasets, svm
    from sklearn.svm import SVR, LinearSVR, SVC, LinearSVC
    from sklearn.kernel_approximation import Nystroem
    from sklearn.tree import DecisionTreeRegressor
    from sklearn.kernel_ridge import KernelRidge
    from sklearn import linear_model
    from sklearn.neural_network import MLPRegressor
    from sklearn.datasets import make_friedman2
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.datasets import make_regression
    from sklearn import ensemble
    from sklearn.utils import shuffle
    from sklearn.metrics import mean_squared_error
    from sklearn.linear_model import LinearRegression
    from sklearn.ensemble import AdaBoostRegressor
    from sklearn.datasets import make_regression
    from sklearn.preprocessing import StandardScaler

    """
 simple_regr_models = [
    GaussianProcessRegressor(),
    RandomForestRegressor(),
    GradientBoostingRegressor(),
    AdaBoostRegressor(),
    SVR(),
    Lasso(),
    LinearRegression(),
    KernelRidge(),
    MLPRegressor(),
    DecisionTreeRegressor(),
    KernelRidge(alpha=0.5,kernel='rbf', gamma=0.1),
    MLPRegressor(activation='logistic',max_iter=1000),
    GaussianProcessRegressor(random_state=0, alpha=1e-6),
    RandomForestRegressor(max_depth=2, random_state=0,n_estimators=100),
    GradientBoostingRegressor(n_estimators=500, max_depth=4, min_samples_split=2,learning_rate=0.01),
    AdaBoostRegressor(random_state=0, n_estimators=100),
    linear_model.Lasso(alpha=7,max_iter=200000,tol=0.001),
 ]
 """

    # Parameters for GradientBoostingRegressor
    # ========================================
    #        "n_estimators": [100],
    #        "loss": ["ls", "lad", "huber", "quantile"],
    #        "learning_rate": [1e-3, 1e-2, 1e-1, 0.5, 1.0],
    #        "max_depth": range(1, 11),
    #        "min_samples_split": range(2, 21),
    #        "min_samples_leaf": range(1, 21),
    #        "subsample": np.arange(0.05, 1.01, 0.05),
    #        "max_features": np.arange(0.05, 1.01, 0.05),
    #        "alpha": [0.75, 0.8, 0.85, 0.9, 0.95, 0.99],

    # clf = GradientBoostingRegressor(loss='quantile', alpha=alpha,
    #                                n_estimators=250, max_depth=3,
    #                               learning_rate=.1, min_samples_leaf=9,
    #                               min_samples_split=9)

    # param_grid = dict(n_estimators=n_estimators,loss=loss,learning_rate=learning_rate,max_depth=max_depth,min_samples_split=min_samples_split,min_samples_leaf=min_samples_leaf,max_features=max_features,alpha=alpha)

    # Set lower and upper quantile
    # 0.5StanDev
    LOWER_ALPHA = 0.33
    MID_ALPHA = 0.50
    UPPER_ALPHA = 0.64
    # StanDev
    LOWER_ALPHA = 0.16
    MID_ALPHA = 0.50
    UPPER_ALPHA = 0.84

    """
 pipe = Pipeline(
     [
       ("stdscal", StandardScaler()),
       ("vart", VarianceThreshold(1e-4)),
       ("est", lower_model),
     ]
 )
 lower_model = pipe
 """

    # SEARCH PARAMETERS
    scoring = "neg_mean_absolute_error"
    cv = 2
    n_jobs = 24
    # n_iter = how many parameter combination to try in the search
    n_iter = 100
    random_state = 508842607

    scaler = StandardScaler().fit(X_train)
    scaler.transform(X_train)
    scaler.transform(X_test)

    # LOWER Model

    # lower_model.fit(X_train, y_train)
    objective = "quantile"
    alpha = LOWER_ALPHA
    print("Prima di lgbm for LOWER model")
    lower_model = get_lgbm(
        X_train,
        X_test,
        y_train,
        y_test,
        cv,
        n_jobs,
        scoring,
        n_iter,
        objective,
        alpha,
        random_state,
    )
    print("Dopo lgbm for LOWER model")
    name = str(prop) + str("_lower")
    filename = str("pickle2-") + str(name) + str(".pk")
    pickle.dump(lower_model, open(filename, "wb"))
    # filename=str('joblib2-')+str(name)+str('.pkl')
    # joblib.dump(lower_model, filename)

    # MID Model

    """
 pipe = Pipeline(
     [
       ("stdscal", StandardScaler()),
       ("vart", VarianceThreshold(1e-4)),
       ("est", mid_model),
     ]
 )
 mid_model = pipe
 """
    scaler = StandardScaler().fit(X_train)
    scaler.transform(X_train)
    scaler.transform(X_test)

    # mid_model.fit(X_train, y_train)
    objective = "regression"
    alpha = 0.9
    print("Prima di lgbm for MID model")
    mid_model = get_lgbm(
        X_train,
        X_test,
        y_train,
        y_test,
        cv,
        n_jobs,
        scoring,
        n_iter,
        objective,
        alpha,
        random_state,
    )
    print("Dopo lgbm for MID model")
    name = str(prop) + str("_mid")
    filename = str("pickle2-") + str(name) + str(".pk")
    pickle.dump(mid_model, open(filename, "wb"))
    # filename=str('joblib2-')+str(name)+str('.pkl')
    # joblib.dump(mid_model, filename)

    # importances = mid_model.feature_importances_
    # indices = np.argsort(importances)[::-1]
    ## Print the feature ranking
    # print("Feature ranking:")
    # for f in range(X_train.shape[1]):
    # print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

    # UPPER Model
    """
 pipe = Pipeline(
     [
       ("stdscal", StandardScaler()),
       ("vart", VarianceThreshold(1e-4)),
       ("est", upper_model),
     ]
 )
 upper_model = pipe
 """
    scaler = StandardScaler().fit(X_train)
    scaler.transform(X_train)
    scaler.transform(X_test)

    # upper_model.fit(X_train, y_train)
    objective = "quantile"
    alpha = UPPER_ALPHA
    print("Prima di lgbm for UPPER model")
    upper_model = get_lgbm(
        X_train,
        X_test,
        y_train,
        y_test,
        cv,
        n_jobs,
        scoring,
        n_iter,
        objective,
        alpha,
        random_state,
    )
    print("Dopo lgbm for UPPER model")
    name = str(prop) + str("_upper")
    filename = str("pickle2-") + str(name) + str(".pk")
    pickle.dump(upper_model, open(filename, "wb"))
    # filename=str('joblib2-')+str(name)+str('.pkl')
    # joblib.dump(upper_model, filename)

    lower = lower_model.predict(X_test)
    mid = mid_model.predict(X_test)
    upper = upper_model.predict(X_test)
    actual = y_test

    print("Model       mae     rmse")
    reg_sc = regr_scores(y_test, lower)
    print("Lower:", round(reg_sc["mae"], 3), round(reg_sc["rmse"], 3))
    reg_sc = regr_scores(y_test, mid)
    print("Mid:", round(reg_sc["mae"], 3), round(reg_sc["rmse"], 3))
    reg_sc = regr_scores(y_test, upper)
    print("Upper:", round(reg_sc["mae"], 3), round(reg_sc["rmse"], 3))

    """
 Calculate the absolute error associated with prediction intervals
 """
    # in_bounds = actual.between(left=lower, right=upper)

    fout1 = open("Intervals.dat", "w")
    fout2 = open("Intervals1.dat", "w")
    line = "#    Jid      Observed       pred_Lower       pred_Mid        pred_Upper\n"
    fout1.write(line)
    line = "#    Jid      Observed       pred_Lower    AbsErr(Lower)     pred_Mid    AbsErr(Mid)     pred_Upper     AbsErr(Upper)    AbsErrInterval    Pred_inBounds\n"
    fout2.write(line)
    sum = 0.0
    for ii in range(len(actual)):
        true = float(actual[ii])
        pred = float(mid[ii])
        llow = float(lower[ii])
        uupper = float(upper[ii])
        if (true >= llow) and (true <= uupper):
            in_bounds = "True"
        else:
            in_bounds = "False"
        absolute_error_lower = abs(lower[ii] - actual[ii])
        absolute_error_mid = abs(mid[ii] - actual[ii])
        absolute_error_upper = abs(upper[ii] - actual[ii])
        absolute_error_interval = (
            absolute_error_lower + absolute_error_upper
        ) / 2.0

        line = (
            str(ii)
            + " "
            + jid[ii]
            + " "
            + str(actual[ii])
            + " "
            + str(lower[ii])
            + " "
            + str(mid[ii])
            + " "
            + str(upper[ii])
            + "\n"
        )
        sum = sum + float(absolute_error_interval)
        line2 = (
            str(ii)
            + " "
            + jid[ii]
            + " "
            + str(actual[ii])
            + " "
            + str(lower[ii])
            + " "
            + str(absolute_error_lower)
            + " "
            + str(mid[ii])
            + " "
            + str(absolute_error_mid)
            + " "
            + str(upper[ii])
            + " "
            + str(absolute_error_upper)
            + " "
            + str(absolute_error_interval)
            + " "
            + str(in_bounds)
            + "\n"
        )
        fout1.write(line)
        fout2.write(line2)
    print("")
    print(
        "Normalized sum of absolute_error_interval= ", sum / float(len(actual))
    )


# epsx
run(prop="formation_energy_peratom")
