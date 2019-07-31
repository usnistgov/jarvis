from jarvis.sklearn.build_models import isfloat, jdata, regr_scores, plot_learning_curve
import numpy as np
import pytest

def test_isfloat():
    assert isfloat(5.5) == True

def test_jdata():
    X, Y, jid = jdata(prop="form_enp")
    assert len(X) == 24759


def test_reg_score():
    x = np.array([1, 2, 3, 4])
    y = np.array([2, 4, 6, 8])
    info = regr_scores(x, y)
    assert info["mae"] == 2.5
