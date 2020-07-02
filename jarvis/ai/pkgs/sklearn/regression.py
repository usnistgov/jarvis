"""
Simple ML models for regression.

Designed for educational purposes only.
"""
from collections import defaultdict
from jarvis.ai.pkgs.utils import regr_scores
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.ensemble import (
    RandomForestRegressor,
    GradientBoostingRegressor,
    AdaBoostRegressor,
)
from sklearn.svm import SVR
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.kernel_ridge import KernelRidge
from sklearn.neural_network import MLPRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt

# Note that these models are with default parameters
# without any hyper-parameter otimizations
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
]


def regression(X=[], Y=[], plot=False, models=simple_regr_models,
               preprocess=True, test_size=0.1):
    """
    Provide model as models to get accuracy.

    Args:
        X: input features

        Y: Target data

        models : collection array of models

        plot: whether to make a parity plot with ML models

        preprocess: whether to apply standard preprocessing techniques
    """
    X_train, X_test, y_train, y_test = train_test_split(
        X, Y, random_state=1, test_size=test_size
    )
    info = defaultdict()
    for i in models:
        pipe = Pipeline(
            [
                ("stdscal", StandardScaler()),
                ("vart", VarianceThreshold(1e-4)),
                ("est", i),
            ]
        )
        if preprocess:
            model = pipe
        else:
            model = i
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        reg_sc = regr_scores(y_test, pred)
        if plot:
            plt.plot(
                reg_sc["pred"], reg_sc["test"],
                ".", label=str(type(i).__name__)[0:4]
            )
        print(type(i).__name__, round(reg_sc["mae"],
              3), round(reg_sc["rmse"], 3))
        info[type(i).__name__] = {}
        info[type(i).__name__]["mae"] = reg_sc["mae"]
        info[type(i).__name__]["rmse"] = reg_sc["rmse"]

    if plot:
        plt.legend()
        plt.xlabel("DFT")
        plt.ylabel("ML")
        plt.show()

    return info


"""
if __name__ == "__main__":
    from jarvis.ai.pkgs.utils import get_ml_data
    X, Y, jid  = get_ml_data(dataset =
                             'cfid_3d', ml_property='exfoliation_energy')
    info = regression(X, Y, models = simple_regr_models)
    print ('info',info)
    assert info['GradientBoostingRegressor']['mae']<50.0
"""
