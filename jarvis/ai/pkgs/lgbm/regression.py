"""Modules for LightGBM regression."""

from sklearn.model_selection import (
    train_test_split,
    RandomizedSearchCV,
)
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline
import lightgbm as lgb
from jarvis.ai.pkgs.utils import regr_scores
from collections import defaultdict
import numpy as np
import pickle
import joblib
import matplotlib.pyplot as plt
import scipy as sp


def regression(
    X=[],
    Y=[],
    jid=[],
    test_size=0.1,
    plot=False,
    preprocess=True,
    feature_importance=True,
    save_model=False,
    feat_names=[],
    model_name="my_model",
    config={},
):
    """Get generic regression model."""
    lgbm = lgb.LGBMRegressor(
        n_estimators=config["n_estimators"],
        learning_rate=config["learning_rate"],
        num_leaves=config["num_leaves"],
    )
    info = defaultdict()

    X_train, X_test, y_train, y_test, jid_train, jid_test = train_test_split(
        X, Y, jid, random_state=1, test_size=test_size
    )
    pipe = Pipeline(
        [
            ("stdscal", StandardScaler()),
            ("vart", VarianceThreshold(1e-4)),
            ("est", lgbm),
        ]
    )
    if preprocess:
        model = pipe
    else:
        model = lgbm
    model.fit(X_train, y_train)
    pred = model.predict(X_test)
    reg_sc = regr_scores(y_test, pred)
    info["reg_scores"] = reg_sc

    if feature_importance:
        imp_data = []
        info["imp_data"] = imp_data
        if not preprocess:
            feat_imp = model.feature_importances_
            feat_imp = 100 * np.array(
                [float(i) / float(np.sum(feat_imp)) for i in feat_imp]
            )
            for f in range(len(feat_imp)):
                imp_data.append([feat_imp[f], feat_names[f]])
        else:
            feat_imp = model.named_steps["est"].feature_importances_
            feat_imp = 100 * np.array(
                [float(i) / float(np.sum(feat_imp)) for i in feat_imp]
            )
            keep_indices = model.named_steps["vart"].get_support(indices=True)
            indices = np.argsort(feat_imp)[::-1]
            new_feat_imp = feat_imp[indices]
            new_indices = keep_indices[indices]
            for f in range(len(new_feat_imp)):
                imp_data.append([new_feat_imp[f], feat_names[new_indices[f]]])

    print(
        model_name,
        round(reg_sc["mae"], 3),
        round(reg_sc["rmse"], 3),
        round(reg_sc["r2"], 3),
    )
    if plot:
        plt.plot(
            reg_sc["pred"],
            reg_sc["test"],
            ".",
            label=str(type(model).__name__)[0:4],
        )
        plt.legend()
        plt.xlabel("DFT")
        plt.ylabel("ML")

    if save_model:
        pk = str(model_name) + str(".pk")
        jb = str(model_name) + str(".jb")
        # js = str(model_name )+str('.js')
        pickle.dump(model, open(pk, "wb"))
        joblib.dump(model, jb)
        # TODO: implemet something like sklearn-json
        # json.dump(model.get_params(), open(js, "w"))
    return info


default_param_dist = {
    # 'boosting_type': [ 'dart'],
    # 'boosting_type': ['gbdt', 'dart', 'rf'],
    # 'num_leaves': sp.stats.randint(2, 1001),
    # 'subsample_for_bin': sp.stats.randint(10, 1001),
    # 'min_split_gain': sp.stats.uniform(0, 5.0),
    # 'min_child_weight': sp.stats.uniform(1e-6, 1e-2),
    # 'reg_alpha': sp.stats.uniform(0, 1e-2),
    # 'reg_lambda': sp.stats.uniform(0, 1e-2),
    # 'tree_learner': ['data', 'feature', 'serial', 'voting' ],
    # 'application': ['regression_l1', 'regression_l2', 'regression'],
    # 'bagging_freq': sp.stats.randint(1, 11),
    # 'bagging_fraction': sp.stats.uniform(.1, 0.9),
    # 'feature_fraction': sp.stats.uniform(.1, 0.9),
    # 'learning_rate': sp.stats.uniform(1e-3, 0.9),
    # 'est__num_leaves': [2,8,16],
    # 'est__min_data_in_leaf': [1,2,4],
    # 'est__learning_rate': [0.005,0.01,0.1],
    # 'est__max_depth': [1,3,5], #sp.stats.randint(1, 501),
    # 'est__n_estimators': [num_iteration,2*num_iteration,5*num_iteration],
    # sp.stats.randint(100, 20001),
    # 'gpu_use_dp': [True, False],
    "est__min_data_in_leaf": sp.stats.randint(5, 20),
    "est__n_estimators": sp.stats.randint(500, 2000),
    "est__num_leaves": sp.stats.randint(100, 500),
    "est__max_depth": sp.stats.randint(8, 20),
    "est__learning_rate": sp.stats.uniform(5e-3, 0.5),
}


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
    param_dist=default_param_dist,
):
    """
    Train a lightgbm model.

    Args:

        train_x: samples used for trainiing

        val_x: validation set

        train_y: train targets

        val_y: validation targets

        cv: # of cross-validations

        n_jobs: for making the job parallel

        scoring: scoring function to use such as MAE

    Returns:
           Best estimator.
    """
    # Get converged boosting iterations with high learning rate,
    # MAE as the convergence crietria
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
    # Increase n_iter for production runs.
    rscv = RandomizedSearchCV(
        estimator=pipe,
        param_distributions=param_dist,
        cv=cv,
        scoring=scoring,
        n_iter=n_iter,
        n_jobs=n_jobs,
        verbose=-1,
        random_state=random_state,
        refit=True,
    )
    model = rscv.fit(train_x, train_y)
    print("Best Score: ", model.best_score_)
    print("Best Params: ", model.best_params_)
    print("Best Estimator: ", model.best_estimator_)
    # return model.best_estimator_
    return model


def parameters_dict():
    """Give example optimized parameters."""
    parameters = {
        "optb88vdw_bandgap": {
            "n_estimators": 324,
            "learning_rate": 0.06414333047469417,
            "num_leaves": 31,
        },
        "mbj_bandgap": {
            "n_estimators": 210,
            "learning_rate": 0.04727272041771037,
            "num_leaves": 121,
        },
        "epsx": {
            "n_estimators": 139,
            "learning_rate": 0.10098329400041395,
            "num_leaves": 527,
        },
        "epsy": {
            "n_estimators": 161,
            "learning_rate": 0.264679564828344,
            "num_leaves": 29,
        },
        "epsz": {
            "n_estimators": 161,
            "learning_rate": 0.264679564828344,
            "num_leaves": 29,
        },
        "mepsx": {
            "n_estimators": 75,
            "learning_rate": 0.05374708509141705,
            "num_leaves": 242,
        },
        "mepsy": {
            "n_estimators": 120,
            "learning_rate": 0.12048289662270327,
            "num_leaves": 398,
        },
        "mepsz": {
            "n_estimators": 89,
            "learning_rate": 0.09718152788954888,
            "num_leaves": 938,
        },
        "encut": {
            "n_estimators": 376,
            "learning_rate": 0.08982089572506267,
            "num_leaves": 762,
        },
        "kpoint_length_unit": {
            "n_estimators": 236,
            "learning_rate": 0.03234907667844313,
            "num_leaves": 794,
        },
        "bulk_modulus_kv": {
            "n_estimators": 380,
            "learning_rate": 0.08621497083536021,
            "num_leaves": 202,
        },
        "shear_modulus_gv": {
            "n_estimators": 284,
            "learning_rate": 0.017555838240950795,
            "num_leaves": 579,
        },
        "formation_energy_peratom": {
            "n_estimators": 1170,
            "learning_rate": 0.15375236057119931,
            "num_leaves": 273,
        },
        "exfoliation_energy": {
            "n_estimators": 47,
            "learning_rate": 0.0734095239247927,
            "num_leaves": 326,
        },
        "max_ir_mode": {
            "n_estimators": 500,
            "learning_rate": 0.0734095239247927,
            "num_leaves": 100,
        },
    }
    return parameters


"""
if __name__ == "__main__":
    from jarvis.ai.pkgs.utils import get_ml_data
    from jarvis.ai.descriptors.cfid import feat_names
    property = "exfoliation_energy"
    property='formation_energy_peratom'
    params = parameters_dict()[property]
    print(params)
    X, Y, jid = get_ml_data(dataset="cfid_3d", ml_property=property)
    names = feat_names()
    info = regression(X=X, Y=Y, jid = jid, config=params, feat_names=names)
    print(info)
"""
