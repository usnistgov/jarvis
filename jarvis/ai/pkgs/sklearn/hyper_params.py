"""
Set of ranges for hyperparameters.

# Modified from  https://github.com/EpistasisLab/tpot
"""

import numpy as np


def classification_regression_params():
    """Set of hyperparameters."""
    regressor_config_dict = {
        "sklearn.linear_model.ElasticNetCV": {
            "l1_ratio": np.arange(0.0, 1.01, 0.05),
            "tol": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
        },
        "sklearn.ensemble.ExtraTreesRegressor": {
            "n_estimators": [100],
            "max_features": np.arange(0.05, 1.01, 0.05),
            "min_samples_split": range(2, 21),
            "min_samples_leaf": range(1, 21),
            "bootstrap": [True, False],
        },
        "sklearn.ensemble.GradientBoostingRegressor": {
            "n_estimators": [100],
            "loss": ["ls", "lad", "huber", "quantile"],
            "learning_rate": [1e-3, 1e-2, 1e-1, 0.5, 1.0],
            "max_depth": range(1, 11),
            "min_samples_split": range(2, 21),
            "min_samples_leaf": range(1, 21),
            "subsample": np.arange(0.05, 1.01, 0.05),
            "max_features": np.arange(0.05, 1.01, 0.05),
            "alpha": [0.75, 0.8, 0.85, 0.9, 0.95, 0.99],
        },
        "sklearn.ensemble.AdaBoostRegressor": {
            "n_estimators": [100],
            "learning_rate": [1e-3, 1e-2, 1e-1, 0.5, 1.0],
            "loss": ["linear", "square", "exponential"],
        },
        "sklearn.tree.DecisionTreeRegressor": {
            "max_depth": range(1, 11),
            "min_samples_split": range(2, 21),
            "min_samples_leaf": range(1, 21),
        },
        "sklearn.neighbors.KNeighborsRegressor": {
            "n_neighbors": range(1, 101),
            "weights": ["uniform", "distance"],
            "p": [1, 2],
        },
        "sklearn.linear_model.LassoLarsCV": {"normalize": [True, False]},
        "sklearn.svm.LinearSVR": {
            "loss": ["epsilon_insensitive", "squared_epsilon_insensitive"],
            "dual": [True, False],
            "tol": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
            "C": [1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1.0,
                  5.0, 10.0, 15.0, 20.0, 25.0],
            "epsilon": [1e-4, 1e-3, 1e-2, 1e-1, 1.0],
        },
        "sklearn.ensemble.RandomForestRegressor": {
            "n_estimators": [100],
            "max_features": np.arange(0.05, 1.01, 0.05),
            "min_samples_split": range(2, 21),
            "min_samples_leaf": range(1, 21),
            "bootstrap": [True, False],
        },
        "sklearn.linear_model.RidgeCV": {},
        "xgboost.XGBRegressor": {
            "n_estimators": [100],
            "max_depth": range(1, 11),
            "learning_rate": [1e-3, 1e-2, 1e-1, 0.5, 1.0],
            "subsample": np.arange(0.05, 1.01, 0.05),
            "min_child_weight": range(1, 21),
            "nthread": [1],
            "objective": ["reg:squarederror"],
        },
        # Preprocesssors
        "sklearn.preprocessing.Binarizer":
        {"threshold": np.arange(0.0, 1.01, 0.05)},
        "sklearn.decomposition.FastICA": {"tol": np.arange(0.0, 1.01, 0.05)},
        "sklearn.cluster.FeatureAgglomeration": {
            "linkage": ["ward", "complete", "average"],
            "affinity": ["euclidean", "l1", "l2", "manhattan", "cosine"],
        },
        "sklearn.preprocessing.MaxAbsScaler": {},
        "sklearn.preprocessing.MinMaxScaler": {},
        "sklearn.preprocessing.Normalizer": {"norm": ["l1", "l2", "max"]},
        "sklearn.kernel_approximation.Nystroem": {
            "kernel": [
                "rbf",
                "cosine",
                "chi2",
                "laplacian",
                "polynomial",
                "poly",
                "linear",
                "additive_chi2",
                "sigmoid",
            ],
            "gamma": np.arange(0.0, 1.01, 0.05),
            "n_components": range(1, 11),
        },
        "sklearn.decomposition.PCA": {
            "svd_solver": ["randomized"],
            "iterated_power": range(1, 11),
        },
        "sklearn.preprocessing.PolynomialFeatures": {
            "degree": [2],
            "include_bias": [False],
            "interaction_only": [False],
        },
        "sklearn.kernel_approximation.RBFSampler": {
            "gamma": np.arange(0.0, 1.01, 0.05)
        },
        "sklearn.preprocessing.RobustScaler": {},
        "sklearn.preprocessing.StandardScaler": {},
        "tpot.builtins.ZeroCount": {},
        "tpot.builtins.OneHotEncoder": {
            "minimum_fraction": [0.05, 0.1, 0.15, 0.2, 0.25],
            "sparse": [False],
            "threshold": [10],
        },
        # Selectors
        "sklearn.feature_selection.SelectFwe": {
            "alpha": np.arange(0, 0.05, 0.001),
            "score_func": {"sklearn.feature_selection.f_regression": None},
        },
        "sklearn.feature_selection.SelectPercentile": {
            "percentile": range(1, 100),
            "score_func": {"sklearn.feature_selection.f_regression": None},
        },
        "sklearn.feature_selection.VarianceThreshold": {
            "threshold": [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2]
        },
        "sklearn.feature_selection.SelectFromModel": {
            "threshold": np.arange(0, 1.01, 0.05),
            "estimator": {
                "sklearn.ensemble.ExtraTreesRegressor": {
                    "n_estimators": [100],
                    "max_features": np.arange(0.05, 1.01, 0.05),
                }
            },
        },
    }

    classifier_config_dict = {
        # Classifiers
        "sklearn.naive_bayes.GaussianNB": {},
        "sklearn.naive_bayes.BernoulliNB": {
            "alpha": [1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0],
            "fit_prior": [True, False],
        },
        "sklearn.naive_bayes.MultinomialNB": {
            "alpha": [1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0],
            "fit_prior": [True, False],
        },
        "sklearn.tree.DecisionTreeClassifier": {
            "criterion": ["gini", "entropy"],
            "max_depth": range(1, 11),
            "min_samples_split": range(2, 21),
            "min_samples_leaf": range(1, 21),
        },
        "sklearn.ensemble.ExtraTreesClassifier": {
            "n_estimators": [100],
            "criterion": ["gini", "entropy"],
            "max_features": np.arange(0.05, 1.01, 0.05),
            "min_samples_split": range(2, 21),
            "min_samples_leaf": range(1, 21),
            "bootstrap": [True, False],
        },
        "sklearn.ensemble.RandomForestClassifier": {
            "n_estimators": [100],
            "criterion": ["gini", "entropy"],
            "max_features": np.arange(0.05, 1.01, 0.05),
            "min_samples_split": range(2, 21),
            "min_samples_leaf": range(1, 21),
            "bootstrap": [True, False],
        },
        "sklearn.ensemble.GradientBoostingClassifier": {
            "n_estimators": [100],
            "learning_rate": [1e-3, 1e-2, 1e-1, 0.5, 1.0],
            "max_depth": range(1, 11),
            "min_samples_split": range(2, 21),
            "min_samples_leaf": range(1, 21),
            "subsample": np.arange(0.05, 1.01, 0.05),
            "max_features": np.arange(0.05, 1.01, 0.05),
        },
        "sklearn.neighbors.KNeighborsClassifier": {
            "n_neighbors": range(1, 101),
            "weights": ["uniform", "distance"],
            "p": [1, 2],
        },
        "sklearn.svm.LinearSVC": {
            "penalty": ["l1", "l2"],
            "loss": ["hinge", "squared_hinge"],
            "dual": [True, False],
            "tol": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
            "C": [1e-4, 1e-3, 1e-2, 1e-1,
                  0.5, 1.0, 5.0, 10.0, 15.0, 20.0, 25.0],
        },
        "sklearn.linear_model.LogisticRegression": {
            "penalty": ["l1", "l2"],
            "C": [1e-4, 1e-3, 1e-2, 1e-1, 0.5,
                  1.0, 5.0, 10.0, 15.0, 20.0, 25.0],
            "dual": [True, False],
        },
        "xgboost.XGBClassifier": {
            "n_estimators": [100],
            "max_depth": range(1, 11),
            "learning_rate": [1e-3, 1e-2, 1e-1, 0.5, 1.0],
            "subsample": np.arange(0.05, 1.01, 0.05),
            "min_child_weight": range(1, 21),
            "nthread": [1],
        },
        # Preprocesssors
        "sklearn.preprocessing.Binarizer":
        {"threshold": np.arange(0.0, 1.01, 0.05)},
        "sklearn.decomposition.FastICA": {"tol": np.arange(0.0, 1.01, 0.05)},
        "sklearn.cluster.FeatureAgglomeration": {
            "linkage": ["ward", "complete", "average"],
            "affinity": ["euclidean", "l1", "l2", "manhattan", "cosine"],
        },
        "sklearn.preprocessing.MaxAbsScaler": {},
        "sklearn.preprocessing.MinMaxScaler": {},
        "sklearn.preprocessing.Normalizer": {"norm": ["l1", "l2", "max"]},
        "sklearn.kernel_approximation.Nystroem": {
            "kernel": [
                "rbf",
                "cosine",
                "chi2",
                "laplacian",
                "polynomial",
                "poly",
                "linear",
                "additive_chi2",
                "sigmoid",
            ],
            "gamma": np.arange(0.0, 1.01, 0.05),
            "n_components": range(1, 11),
        },
        "sklearn.decomposition.PCA": {
            "svd_solver": ["randomized"],
            "iterated_power": range(1, 11),
        },
        "sklearn.preprocessing.PolynomialFeatures": {
            "degree": [2],
            "include_bias": [False],
            "interaction_only": [False],
        },
        "sklearn.kernel_approximation.RBFSampler": {
            "gamma": np.arange(0.0, 1.01, 0.05)
        },
        "sklearn.preprocessing.RobustScaler": {},
        "sklearn.preprocessing.StandardScaler": {},
        "tpot.builtins.ZeroCount": {},
        "tpot.builtins.OneHotEncoder": {
            "minimum_fraction": [0.05, 0.1, 0.15, 0.2, 0.25],
            "sparse": [False],
            "threshold": [10],
        },
        # Selectors
        "sklearn.feature_selection.SelectFwe": {
            "alpha": np.arange(0, 0.05, 0.001),
            "score_func": {"sklearn.feature_selection.f_classif": None},
        },
        "sklearn.feature_selection.SelectPercentile": {
            "percentile": range(1, 100),
            "score_func": {"sklearn.feature_selection.f_classif": None},
        },
        "sklearn.feature_selection.VarianceThreshold": {
            "threshold": [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2]
        },
        "sklearn.feature_selection.RFE": {
            "step": np.arange(0.05, 1.01, 0.05),
            "estimator": {
                "sklearn.ensemble.ExtraTreesClassifier": {
                    "n_estimators": [100],
                    "criterion": ["gini", "entropy"],
                    "max_features": np.arange(0.05, 1.01, 0.05),
                }
            },
        },
        "sklearn.feature_selection.SelectFromModel": {
            "threshold": np.arange(0, 1.01, 0.05),
            "estimator": {
                "sklearn.ensemble.ExtraTreesClassifier": {
                    "n_estimators": [100],
                    "criterion": ["gini", "entropy"],
                    "max_features": np.arange(0.05, 1.01, 0.05),
                }
            },
        },
    }
    return regressor_config_dict, classifier_config_dict
