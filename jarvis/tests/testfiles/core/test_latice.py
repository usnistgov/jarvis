from jarvis.core.lattice import Lattice
import numpy as np


def test_lat():
    box = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    lat = Lattice(box)
    frac_coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
    cart_coords = [[0, 0, 0], [5, 5, 5]]
    lll = lat._calculate_lll()
    # print ('lll',lll[0][0][0])
    lll_red = lat.get_lll_reduced_lattice()._lat
    # print("lll_educed", lat.get_lll_reduced_lattice()._lat[0][0])
    assert (
        lat.lat_lengths(),
        lat.lat_angles(),
        round(lat.inv_lattice()[0][0], 2),
        [round(i, 2) for i in lat.lat_angles(radians=True)],
        lat.cart_coords(frac_coords)[1][1],
        lat.frac_coords(cart_coords)[1][1],
        lat.volume,
        lat.parameters,
        lll[0][0][0],
        lll_red[0][0],
    ) == (
        [10.0, 10.0, 10.0],
        [90.0, 90.0, 90.0],
        0.1,
        [1.57, 1.57, 1.57],
        5.0,
        0.5,
        1000.0,
        [10.0, 10.0, 10.0, 90.0, 90.0, 90.0],
        10.0,
        10.0,
    )


# test_lat()
