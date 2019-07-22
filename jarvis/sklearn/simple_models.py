"""
Simple ML models for classifcation and regression, designed for educational purposes only
__author__: = Kamal Choudhary
"""

import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.ensemble import (
    RandomForestRegressor,
    GradientBoostingRegressor,
    AdaBoostRegressor,
)
from sklearn.svm import SVR
from sklearn.linear_model import Lasso, LinearRegression, LogisticRegression
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
from pymatgen.core.structure import Structure
from jarvis.sklearn.build_models import *
from jarvis.sklearn.build_models import jdata
from jarvis.db.static.explore_db import get_3d_dataset, get_2d_dataset, get_ml_dataset

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

simple_class_models = [
    DecisionTreeClassifier(),
    MLPClassifier(),
    GradientBoostingClassifier(),
    KNeighborsClassifier(),
    GaussianProcessClassifier(),
    RandomForestClassifier(),
    AdaBoostClassifier(),
    SVC(),
]


def read_vasp_structure(filename=""):
    """
   Read VASP POSCAR, .cif file
   Args:
       filename:string
   Returns:
       s: pymatgen.structure object
   """

    s = Structure.from_file(filename)
    return s


def simple_regression(
    X=[], Y=[], plot=False, simple_models=simple_regr_models, preprocess=True
):
    """
  Quickly train simple regression models without hyperparameter-optimization
  Args:
      X: input features
      Y: Target data
      plot: whether to make a parity plot with ML models
      preprocess: whether to apply standard preprocessing techniques
  Returns:
  """

    X_train, X_test, y_train, y_test = train_test_split(
        X, Y, random_state=1, test_size=0.1
    )
    for i in simple_models:
        pipe = Pipeline(
            [
                ("stdscal", StandardScaler()),
                ("vart", VarianceThreshold(1e-4)),
                ("est", i),
            ]
        )
        if preprocess == True:
            model = pipe
        else:
            model = i
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        reg_sc = regr_scores(y_test, pred)
        if plot == True:
            plt.plot(
                reg_sc["pred"], reg_sc["test"], ".", label=str(type(i).__name__)[0:4]
            )
        print(type(i).__name__, round(reg_sc["mae"], 3), round(reg_sc["rmse"], 3))
    if plot == True:
        plt.legend()
        plt.xlabel("DFT")
        plt.ylabel("ML")


def binary_class_dat(X=[], Y=[], tol=0.1):
    """
    Categorize a continous dataset in 1/0 with a threshold "tol"
    """

    Y1 = []
    for i, j in zip(X, Y):
        if j >= tol:
            Y1.append(1)
        else:
            Y1.append(0)
    return X, Y1


def bar_plot(x=[], interval=10):
    """
  Plot bar plot to understand the range of data
  """
    # interval = # abs(max(x)-min(x))/len(x)
    hist, bins = np.histogram(
        x, bins=np.arange(min(x), max(x), interval), density=False
    )
    return bins[:-1], hist


def classify_roc_ml(
    X=[],
    y=[],
    classes=[0, 1, 2],
    names=["High val", "Low val", ""],
    n_plot=1,
    method="",
    preprocess=True,
    plot=False,
):
    """
    Classifcation module for ROC curve for upto three classes, can be expanded in more classes as well
    Args:
        X: input feature vectors
        y: target data obtained from binary_class_dat
        classes: dummy classes
        names: name holders for the target data
        method: ML method
        preprocess: whether to apply standard preprocessing techniques
        plt: whether to plot the ROC curve
    """
    if plot == True:
        plt.close()
        plt.rcParams.update({"font.size": 22})
        plt.figure(figsize=(12, 8))

    y = label_binarize(y, classes=classes)
    n_classes = y.shape[1]

    pipe = Pipeline(
        [
            ("stdscal", StandardScaler()),
            ("vart", VarianceThreshold(1e-4)),
            ("est", method),
        ]
    )
    if preprocess == True:
        model = pipe
    else:
        model = method
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.1, random_state=0
    )
    classifier = OneVsRestClassifier(model)
    if hasattr(model, "decision_function"):
        y_score = classifier.fit(X_train, y_train).decision_function(X_test)
    else:
        y_score = classifier.fit(X_train, y_train).predict_proba(X_test)

    lw = 3

    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])
    colors = ["blue", "red", "green"]
    count = 0
    for i, color, name in zip(range(n_classes), colors, names):

        if name != "":
            if count < n_plot:
                count = count + 1
                if plot == True:
                    plt.plot(
                        fpr[i],
                        tpr[i],
                        color=color,
                        lw=lw,
                        label="ROC  {0} (area = {1:0.2f})" "".format(name, roc_auc[i]),
                    )
    if plot == True:
        plt.plot([0, 1], [0, 1], "k--", lw=lw)
        plt.xlim([-0.05, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.legend(loc="lower right")
        plt.show()
    model = classifier.fit(X, y)
    return model, roc_auc


def simple_classification(X=[],Y=[],tol=100,plot=True, simple_models=simple_class_models):
    """
    Quickly train some of the classifcation algorithms available in scikit-learn
    """
    X_class, Y_class = binary_class_dat(X=X, Y=Y, tol=tol)
    for i in simple_models:
        m, r = classify_roc_ml(X=X_class, y=Y_class, method=i, plot=plot)
        print(type(i).__name__, r[0])


def simple_GB_confidence_interval(X=[], Y=[]):
    """
 From https://towardsdatascience.com/how-to-generate-prediction-intervals-with-scikit-learn-and-python-ab3899f992ed

 """
    LOWER_ALPHA = 0.1
    UPPER_ALPHA = 0.9
    # Each model has to be separate
    lower_model = GradientBoostingRegressor(loss="quantile", alpha=LOWER_ALPHA)
    # The mid model will use the default loss
    mid_model = GradientBoostingRegressor(loss="ls")
    upper_model = GradientBoostingRegressor(loss="quantile", alpha=UPPER_ALPHA)

    models = [lower_model, mid_model, upper_model]
    simple_regression(X=X, Y=Y, simple_models=models)


if __name__ == "__main__":
    X, Y, jid = jdata(prop="exfoliation_en")
    #simple_GB_confidence_interval(X, Y)
    #simple_regression(X, Y)
    #X_class, Y_class = binary_class_dat(X=X, Y=Y, tol=100)
    simple_classification(X=X,Y=Y,tol=100)
