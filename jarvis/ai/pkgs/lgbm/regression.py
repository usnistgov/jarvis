from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline
import lightgbm as lgb
from jarvis.ai.pkgs.utils import get_ml_data
from jarvis.ai.pkgs.utils import regr_scores
from jarvis.ai.descriptors.cfid import feat_names
from collections import defaultdict
import numpy as np
import pickle, joblib, json


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
    if preprocess == True:
        model = pipe
    else:
        model = lgbm
    model.fit(X_train, y_train)
    pred = model.predict(X_test)
    reg_sc = regr_scores(y_test, pred)
    info["reg_scores"] = reg_sc

    if feature_importance == True:
        imp_data = []
        info["imp_data"] = imp_data
        if preprocess != True:
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
    if plot == True:
        plt.plot(reg_sc["pred"], reg_sc["test"], ".", label=str(type(i).__name__)[0:4])
        plt.legend()
        plt.xlabel("DFT")
        plt.ylabel("ML")

    if save_model == True:
        pk = str(model_name) + str(".pk")
        jb = str(model_name) + str(".jb")
        # js = str(model_name )+str('.js')
        pickle.dump(model, open(pk, "wb"))
        joblib.dump(model, jb)
        # TODO: implemet something like sklearn-json
        # json.dump(model.get_params(), open(js, "w"))
    return info


def parameters_dict():
    """Example optimized parameters"""

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
            #"num_leaves": 73,
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


#"""
if __name__ == "__main__":
    property = "exfoliation_energy"
    property='formation_energy_peratom'
    params = parameters_dict()[property]
    print(params)
    X, Y, jid = get_ml_data(dataset="cfid_3d", ml_property=property)
    names = feat_names()
    info = regression(X=X, Y=Y, jid = jid, config=params, feat_names=names)
    print(info)
#"""
