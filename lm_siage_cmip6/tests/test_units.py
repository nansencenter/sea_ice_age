import numpy as np

from lm_siage_cmip6.io import concentration_fraction, velocity_km_day


def test_concentration_fraction_masks_invalid_values():
    values = np.ma.array([0.0, 15.0, 100.0, 1.0e20], mask=[0, 0, 0, 1])
    result = concentration_fraction(values)
    np.testing.assert_allclose(result.compressed(), [0.0, 0.15, 1.0])
    assert result.mask.tolist() == [False, False, False, True]


def test_velocity_km_day():
    result = velocity_km_day(np.array([0.0, 1.0, -0.5]))
    np.testing.assert_allclose(result, [0.0, 86.4, -43.2])
