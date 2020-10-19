"""
Simple ML models for classifcation and regression.

Designed for educational purposes only
"""
from collections import defaultdict
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
from sklearn.multiclass import OneVsRestClassifier
from sklearn.ensemble import (
    RandomForestClassifier,
    AdaBoostClassifier,
    GradientBoostingClassifier,
)
from sklearn.svm import SVC
import pickle
import joblib
import matplotlib.pyplot as plt
from jarvis.ai.pkgs.utils import binary_class_dat

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


def classify_roc_ml(
    X=[],
    y=[],
    classes=[0, 1, 2],
    names=["High val", "Low val", ""],
    n_plot=1,
    method="",
    preprocess=True,
    plot=False,
    test_size=0.1,
):
    """
    Classifcation module for ROC curve for upto three classes.

    It can be expanded in more classes as well.
    Args:
        X: input feature vectors

        y: target data obtained from binary_class_dat

        classes: dummy classes

        names: name holders for the target data

        method: ML method

        preprocess: whether to apply standard preprocessing techniques

        plot: whether to plot the ROC curve
    """
    if plot:
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
    if preprocess:
        model = pipe
    else:
        model = method
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=0
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
                if plot:
                    plt.plot(
                        fpr[i],
                        tpr[i],
                        color=color,
                        lw=lw,
                        label="ROC  {0} (area = {1:0.2f})"
                        "".format(name, roc_auc[i]),
                    )
    if plot:
        plt.plot([0, 1], [0, 1], "k--", lw=lw)
        plt.xlim([-0.05, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.legend(loc="lower right")
        plt.show()
    model = classifier.fit(X, y)
    return model, roc_auc


def classification(
    X=[],
    Y=[],
    tol=100,
    plot=False,
    preprocess=True,
    models=simple_class_models,
    model_name="my_model",
    save_model=False,
):
    """Quickly train some of the classifcation algorithms in scikit-learn."""
    X_class, Y_class = binary_class_dat(X=X, Y=Y, tol=tol)
    info = defaultdict()
    for i in models:
        m, r = classify_roc_ml(
            X=X_class, y=Y_class, method=i, preprocess=preprocess, plot=plot
        )
        print(type(i).__name__, r[0])
        info[type(i).__name__] = {}
        info[type(i).__name__]["roc_auc"] = r
        if save_model:
            pk = (
                str(model_name)
                + "_"
                + str(type(i).__name__)
                + "_"
                + str(".pk")
            )
            jb = (
                str(model_name)
                + "_"
                + str(type(i).__name__)
                + "_"
                + str(".jb")
            )
            pickle.dump(m, open(pk, "wb"))
            joblib.dump(m, jb)
    return info


"""
if __name__ == "__main__":
    from jarvis.ai.pkgs.utils import get_ml_data
    X, Y, jid  = get_ml_data(dataset = 'cfid_3d',
                             ml_property='exfoliation_energy')
    X_class, Y_class = binary_class_dat(X=X, Y=Y, tol=100)
    info = classification(X=X,Y=Y,tol=100, save_model=True)
    #print (info)
    #print ()
    #print (info['GradientBoostingClassifier']['roc_auc'][0])
"""
