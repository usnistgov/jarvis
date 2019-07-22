from jarvis.sklearn.get_desc import get_comp_descp
from monty.serialization import loadfn, MontyDecoder, dumpfn
import numpy as np
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error

try:
    import lightgbm.LGBMRegressor as lgb
except:
    print("WARNING!!!!: LightGBM is not installed, using sklearn GBM")
    from sklearn.ensemble import GradientBoostingRegressor as lgb

import matplotlib.pyplot as plt
from jarvis.db.static.explore_db import get_ml_dataset

# import lightgbm as lgb
import matplotlib.pyplot as plt


plt.switch_backend("agg")
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

"""
Steps for building a regression model
Download descriptors and material data from the following link:
https://figshare.com/articles/JARVIS-ML-CFID-descriptors_and_material_properties/6870101 
You can also make descriptors using get_comp_descp
Link to the paper: https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.083801
Websites: https://www.ctcms.nist.gov/jarvisml, https://jarvis.nist.gov
"""


def isfloat(value):
    """
  Check if an argument is float
  
  Args:
      value: string/int/float
  Returns:
         True/False 
  """

    try:
        float(value)
        return True
    except:
        return False
        pass


def jdata(data_file="jarvisml_cfid.json", prop=""):
    """
  Get data and descriptors for a particular properties
 
  Args:
      data_file: JSON containing descriptors and target properties
      It can be better to save JSON rather tha converting materials
      to descriptors all the time
      prop: property to train such as formation energy, bandgap etc.
  Returns:
         Generic X, Y and identifiers 
         The IDs are used to detect which samples could be trained well and
         which couldn't

  """

    d3 = get_ml_dataset()

    X = []
    Y = []
    jid = []
    for ii, i in enumerate(d3):
        y = i[prop]
        if isfloat(y):
            y = float(y)
            x = i["desc"]
            if len(x) == 1557 and any(np.isnan(x) for x in x.flatten()) == False:
                if "eps" in prop:
                    y = np.sqrt(float(y))
                if "mag" in prop:
                    num = get_number_formula_unit(i["strt"])
                    y = float(abs(y)) / float(num)
                X.append(x)
                Y.append(y)
                jid.append(i["jid"])
    print("Prop=", prop, len(X), len(Y))
    X = np.array(X).astype(np.float64)
    Y = np.array(Y).astype(np.float64)
    return X, Y, jid


def plot_learning_curve(
    estimator=None,
    scoring="neg_mean_absolute_error",
    title="",
    X=None,
    y=None,
    ylim=None,
    cv=5,
    # n_jobs=-1, train_sizes=np.linspace(.1, 1.0, 10),fname='fig.png'):
    n_jobs=1,
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
        estimator, X, y, cv=cv, n_jobs=n_jobs, train_sizes=train_sizes, scoring=scoring
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
    plt.plot(train_sizes, train_scores_mean, "o-", color="r", label="Training score")
    plt.plot(
        train_sizes, test_scores_mean, "o-", color="g", label="Cross-validation score"
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
        train_sizes, test_scores_mean, "o-", color="g", label="Cross-validation score"
    )
    plt.legend(loc="best")
    fname = str(fname) + str("_learn_2.png")
    plt.show()
    plt.savefig(fname)
    plt.close()
    return mem


def regr_scores(pred, test):

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


def get_lgbm(train_x, val_x, train_y, val_y, cv, n_jobs, scoring):
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
           Best estomator
    """

    # Get converged boosting iterations with high learning rate, MAE as the convergence crietria

    lgbm = lgb(
        n_estimators=1000,
        learning_rate=0.1,
        max_depth=5,
        num_leaves=100,
        objective="regression",
        # min_data_in_leaf=2,
        n_jobs=-1,
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
        #'est__num_leaves': sp.stats.randint(3, 1000),
        #'est__max_depth': sp.stats.randint(1, 5),
        "est__learning_rate": sp.stats.uniform(1e-3, 0.9)
    }

    lgbm = lgb(
        objective="regression",
        # device='gpu',
        n_estimators=num_iteration,
        n_jobs=n_jobs,
        verbose=-1,
    )
    pipe = Pipeline(
        [
            ("stdscal", StandardScaler()),
            ("vart", VarianceThreshold(1e-4)),
            ("est", lgbm),
        ]
    )

    n_iter = 10
    # Increase n_iter
    rscv = RandomizedSearchCV(
        estimator=pipe,
        param_distributions=param_dist,
        cv=cv,
        scoring=scoring,
        n_iter=n_iter,
        n_jobs=n_jobs,
        verbose=3,
        refit=True,
    )
    rscv = rscv.fit(train_x, train_y)
    return rscv.best_estimator_


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
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    os.chdir(dir_name)

    info = {}
    tmp_time = time.time()

    # STEP-1: Data for a particular model

    x, y, jid = jdata(prop=prop)

    # Toy example with boston data
    # boston=load_boston()
    # x, y = boston['data'], boston['target']

    # STEP-2: Splitting the data
    # 90-10% split for train test
    X_train, X_test, y_train, y_test, jid_train, jid_test = train_test_split(
        x, y, jid, random_state=1, test_size=0.1
    )

    # further split for boosting iteration convergence
    X_train1, X_test1, y_train1, y_test1 = train_test_split(
        X_train, y_train, random_state=1, test_size=0.1
    )
    print("lenx len y", len(x[0]), len(y))

    # STEP-3: GBM model with converged iterations

    model = get_lgbm(X_train1, X_test1, y_train1, y_test1, cv, n_jobs, scoring)

    model.fit(X_train, y_train)
    print(model, model.named_steps["est"].n_estimators)
    print("model", model.__dict__)
    # info['model']=model.__dict__

    # STEP-4: Predict on 10% held data and  accuracy

    pred = model.predict(X_test)
    reg_sc = regr_scores(y_test, pred)

    print("reg_sc", reg_sc)
    info["reg_sc"] = reg_sc
    info["y"] = y
    info["time"] = time.time() - tmp_time
    info["jid_test"] = jid_test

    model.fit(X_train1, y_train1)
    print("Plot feature importances...")

    # STEP-5: Feature importance
    feat_imp = model.named_steps["est"].feature_importances_  # feature_importances_
    feat_imp = np.array([float(i) / float(np.sum(feat_imp)) for i in feat_imp])
    print("feat_imp", name, feat_imp, len(feat_imp))

    # Since variance threshold removed some features
    keep_indices2 = model.named_steps["vart"].get_support(indices=True)
    info["keep_indices_sel"] = keep_indices2

    print("keep_indices2", keep_indices2, len(keep_indices2))

    indices = np.argsort(feat_imp)[::-1]
    xs = feat_imp[indices]
    ys = keep_indices2[indices]

    print("xs,ys", len(xs), len(ys))
    data_imp = []
    f = open("feat_names.json", "r")
    names = json.load(f)
    f.close()

    for f in range(len(xs)):
        print(f, xs[f], names[ys[f]])
        data_imp.append([f, xs[f], names[ys[f]]])
    info["feature_importance"] = data_imp
    # STEP-6: N-fold cross-validation prediction

    scores = cross_val_score(
        model, x, y, cv=cv, n_jobs=n_jobs, scoring="mean_absolute_error"
    )
    print("MAE-Score: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
    info["cross_val_scores"] = scores

    fname = str(name) + str("-GB-CV.png")
    predicted = cross_val_predict(model, x, y, cv=cv, n_jobs=n_jobs)
    fig, ax = plt.subplots()
    ax.set_xlabel("DFT prediction")
    ax.set_ylabel("ML prediction")
    print("cv_predicted")
    cvpred = {}
    cvpred["cross_val_y"] = y
    cvpred["cross_val_predicted"] = predicted
    info["cvpredict_result"] = cvpred

    ax.scatter(y, predicted, edgecolors=(0, 0, 0))
    plt.savefig(fname)
    plt.close()

    # STEP-7: Plot learning curve

    mem_learn = plot_learning_curve(
        estimator=model, scoring=scoring, fname=name, title="", X=x, y=y, n_jobs=n_jobs
    )
    info["learning_curve_data"] = mem_learn

    # STEP-8: Store ML parameters with pickle and joblib

    model.fit(x, y)
    filename = str("pickle2-") + str(name) + str(".pk")
    pickle.dump(model, open(filename, "wb"))
    filename = str("joblib2-") + str(name) + str(".pkl")
    joblib.dump(model, filename)

    # STEP-9: Store data in a json file

    file = str(name) + "_info.json"
    f = open(file, "w")
    f.write(json.dumps(info, cls=MontyEncoder, indent=4))
    f.close()
    os.chdir("../")


if __name__ == "__main__":
    # This may take long time
    # run(version='version_1',scoring='neg_mean_absolute_error',cv=5,n_jobs=1,prop='op_gap',do_cv=False)

    # smaller test fit model
    model = lgb(
        n_estimators=100,
        learning_rate=0.1,
        max_depth=5,
        num_leaves=100,
        objective="regression",
        n_jobs=-1,
        verbose=-1,
    )
    x, y, jid = jdata(prop="form_enp")
    X_train, X_test, y_train, y_test, jid_train, jid_test = train_test_split(
        x, y, jid, random_state=1, test_size=0.1
    )
    len(X_train), len(X_test)

    # Let's take 500 of training set as a quick example
    X = X_train[0:500]
    Y = y_train[0:500]
    model.fit(X, Y)

    info = {}
    tmp_time = time.time()

    # STEP-1: Data for a particular model

    x, y, jid = jdata(prop=prop)

    # Toy example with boston data
    # boston=load_boston()
    # x, y = boston['data'], boston['target']

    # STEP-2: Splitting the data
    # 90-10% split for train test
    X_train, X_test, y_train, y_test, jid_train, jid_test = train_test_split(
        x, y, jid, random_state=1, test_size=0.1
    )

    # further split for boosting iteration convergence
    X_train1, X_test1, y_train1, y_test1 = train_test_split(
        X_train, y_train, random_state=1, test_size=0.1
    )
    print("lenx len y", len(x[0]), len(y))

    # STEP-3: GBM model with converged iterations

    model = get_lgbm(X_train1, X_test1, y_train1, y_test1, cv, n_jobs, scoring)

    model.fit(X_train, y_train)
    print(model, model.named_steps["est"].n_estimators)
    print("model", model.__dict__)
    # info['model']=model.__dict__

    # STEP-4: Predict on 10% held data and  accuracy

    pred = model.predict(X_test)
    reg_sc = regr_scores(y_test, pred)

    print("reg_sc", reg_sc)
    info["reg_sc"] = reg_sc
    info["y"] = y
    info["time"] = time.time() - tmp_time
    info["jid_test"] = jid_test

    model.fit(X_train1, y_train1)
    print("Plot feature importances...")

    # STEP-5: Feature importance
    feat_imp = model.named_steps["est"].feature_importances_  # feature_importances_
    feat_imp = np.array([float(i) / float(np.sum(feat_imp)) for i in feat_imp])
    print("feat_imp", name, feat_imp, len(feat_imp))

    # Since variance threshold removed some features
    keep_indices2 = model.named_steps["vart"].get_support(indices=True)
    info["keep_indices_sel"] = keep_indices2

    print("keep_indices2", keep_indices2, len(keep_indices2))

    indices = np.argsort(feat_imp)[::-1]
    xs = feat_imp[indices]
    ys = keep_indices2[indices]

    print("xs,ys", len(xs), len(ys))
    data_imp = []
    f = open("feat_names.json", "r")
    names = json.load(f)
    f.close()

    for f in range(len(xs)):
        print(f, xs[f], names[ys[f]])
        data_imp.append([f, xs[f], names[ys[f]]])
    info["feature_importance"] = data_imp
    # STEP-6: N-fold cross-validation prediction

    scores = cross_val_score(
        model, x, y, cv=cv, n_jobs=n_jobs, scoring="mean_absolute_error"
    )
    print("MAE-Score: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
    info["cross_val_scores"] = scores

    fname = str(name) + str("-GB-CV.png")
    predicted = cross_val_predict(model, x, y, cv=cv, n_jobs=n_jobs)
    fig, ax = plt.subplots()
    ax.set_xlabel("DFT prediction")
    ax.set_ylabel("ML prediction")
    print("cv_predicted")
    cvpred = {}
    cvpred["cross_val_y"] = y
    cvpred["cross_val_predicted"] = predicted
    info["cvpredict_result"] = cvpred

    ax.scatter(y, predicted, edgecolors=(0, 0, 0))
    plt.savefig(fname)
    plt.close()

    # STEP-7: Plot learning curve

    mem_learn = plot_learning_curve(
        estimator=model, scoring=scoring, fname=name, title="", X=x, y=y, n_jobs=n_jobs
    )
    info["learning_curve_data"] = mem_learn

    # STEP-8: Store ML parameters with pickle and joblib

    model.fit(x, y)
    filename = str("pickle2-") + str(name) + str(".pk")
    pickle.dump(model, open(filename, "wb"))
    filename = str("joblib2-") + str(name) + str(".pkl")
    joblib.dump(model, filename)

    # STEP-9: Store data in a json file

    file = str(name) + "_info.json"
    f = open(file, "w")
    f.write(json.dumps(info, cls=MontyEncoder, indent=4))
    f.close()
    os.chdir("../")


if __name__ == "__main__":
    # This may take long time
    # run(version='version_1',scoring='neg_mean_absolute_error',cv=5,n_jobs=1,prop='op_gap',do_cv=False)

    # smaller test fit model
    model = lgb(
        n_estimators=100,
        learning_rate=0.1,
        max_depth=5,
        num_leaves=100,
        objective="regression",
        n_jobs=-1,
        verbose=-1,
    )
    x, y, jid = jdata(prop="form_enp")
    X_train, X_test, y_train, y_test, jid_train, jid_test = train_test_split(
        x, y, jid, random_state=1, test_size=0.1
    )
    len(X_train), len(X_test)

    # Let's take 500 of training set as a quick example
    X = X_train[0:500]
    Y = y_train[0:500]
    model.fit(X, Y)

    pred = model.predict(X_test)
    reg_sc = regr_scores(y_test, pred)
    print (reg_sc)
