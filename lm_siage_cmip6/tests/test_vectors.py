import numpy as np

from lm_siage_cmip6.vectors import east_north_to_projected


LAEA = "+proj=laea +lat_0=90 +lon_0=0 +datum=WGS84 +units=m +no_defs"


def test_basis_orientation_on_prime_meridian():
    dx, dy = east_north_to_projected(
        np.array([0.0]), np.array([80.0]),
        np.array([1.0]), np.array([0.0]), LAEA,
    )
    assert dx[0] > 0.9
    assert abs(dy[0]) < 1.0e-3

    dx, dy = east_north_to_projected(
        np.array([0.0]), np.array([80.0]),
        np.array([0.0]), np.array([1.0]), LAEA,
    )
    assert abs(dx[0]) < 1.0e-3
    assert dy[0] > 0.9


def test_rotation_at_ninety_degrees_east():
    dx, dy = east_north_to_projected(
        np.array([90.0]), np.array([80.0]),
        np.array([1.0]), np.array([0.0]), LAEA,
    )
    assert abs(dx[0]) < 1.0e-3
    assert dy[0] > 0.9


def test_masked_fill_values_remain_masked():
    eastward = np.ma.array([0.1, 1.0e20], mask=[False, True])
    northward = np.ma.array([0.0, 1.0e20], mask=[False, True])
    dx, dy = east_north_to_projected(
        np.array([0.0, 0.0]), np.array([80.0, 80.0]),
        eastward, northward, LAEA,
    )
    assert not dx.mask[0]
    assert dx.mask[1]
    assert dy.mask[1]
