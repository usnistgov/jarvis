from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline
from jarvis.sklearn.build_models import jdata
import lightgbm as lgb


def regression(X=[], Y=[]):
    lgbm = lgb.LGBMRegressor(
        n_estimators=100,
        learning_rate=0.1,
        max_depth=5,
        num_leaves=100,
        objective="regression",
        n_jobs=-1,
        verbose=-1,
    )

    X_train, X_test, y_train, y_test, jid_train, jid_test = train_test_split(
        X, Y, jid, random_state=1, test_size=0.1
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
    if plot == True:
        plt.plot(reg_sc["pred"], reg_sc["test"], ".", label=str(type(i).__name__)[0:4])
    print(type(i).__name__, round(reg_sc["mae"], 3), round(reg_sc["rmse"], 3))
    if plot == True:
        plt.legend()
        plt.xlabel("DFT")
        plt.ylabel("ML")


if __name__ == "__main__":
    X, Y, jid = jdata(prop="exfoliation_en")
    regression(X=X, Y=Y)
