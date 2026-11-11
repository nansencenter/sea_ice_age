import numpy as np

from lm_siage_cmip6.interpolation import sample_sic_to_elements


def test_model_sic_interpolation_renormalizes_over_ocean_cells():
    concentration = np.ma.array(
        [[80.0, 80.0], [80.0, 80.0]],
        mask=[[False, True], [False, True]],
    )
    land = np.array([[False, True], [False, True]])
    points = np.array([[0.5, 0.25]])

    model, missing = sample_sic_to_elements(
        np.array([0.0, 1.0]), np.array([0.0, 1.0]),
        concentration, land, points, mode="model_linear",
    )
    legacy, _ = sample_sic_to_elements(
        np.array([0.0, 1.0]), np.array([0.0, 1.0]),
        concentration, land, points, mode="legacy_zero_fill",
    )

    np.testing.assert_allclose(model, [80.0])
    np.testing.assert_allclose(legacy, [60.0])
    assert not missing.any()


def test_model_sic_interpolation_reports_uncovered_ocean():
    concentration = np.ma.masked_all((2, 2))
    sic, missing = sample_sic_to_elements(
        np.array([0.0, 1.0]), np.array([0.0, 1.0]),
        concentration, np.zeros((2, 2), dtype=bool),
        np.array([[0.5, 0.5]]), mode="model_linear",
    )
    assert missing.all()
    assert np.isnan(sic).all()


def test_sic_is_zero_only_when_element_centroid_is_land():
    concentration = np.ma.array(
        [[80.0, 0.0], [80.0, 0.0]],
        mask=[[False, True], [False, True]],
    )
    land = np.array([[False, True], [False, True]])
    sic, missing = sample_sic_to_elements(
        np.array([0.0, 1.0]), np.array([0.0, 1.0]),
        concentration, land, np.array([[0.5, 0.9]]), mode="model_linear",
    )
    np.testing.assert_allclose(sic, [0.0])
    assert not missing.any()
