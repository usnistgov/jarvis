from jarvis.db.static.explore_db import get_3d_dataset, get_2d_dataset, get_ml_dataset
import pytest

def test_j2d():
    data_2d = get_2d_dataset()
    assert (len(data_2d)) == 754

def test_j3d():
    data_3d = get_3d_dataset()
    assert (len(data_3d)) == 32486

def test_jml():
    data_ml = get_ml_dataset()
    assert (len(data_ml)) == 24759
