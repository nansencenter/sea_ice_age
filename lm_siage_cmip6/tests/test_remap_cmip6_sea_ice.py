from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).parents[1] / "scripts" / "remap_cmip6_sea_ice.py"
SPEC = spec_from_file_location("remap_cmip6_sea_ice", SCRIPT)
MODULE = module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def completion_case():
    coordinates = np.arange(3, dtype=np.float64) * 25_000.0
    x_target, y_target = np.meshgrid(coordinates, coordinates)
    lon_target = np.zeros((3, 3), dtype=np.float64)
    lat_target = np.full((3, 3), 80.0, dtype=np.float64)
    target_valid = np.ones((3, 3), dtype=bool)
    concentration = np.zeros((3, 3), dtype=np.float64)
    concentration[1, 1] = 20.0
    u_target = np.full((3, 3), np.nan, dtype=np.float64)
    v_target = np.full((3, 3), np.nan, dtype=np.float64)
    u_target[1, 0] = 2.0
    v_target[1, 0] = -1.0
    return (
        u_target,
        v_target,
        concentration,
        x_target,
        y_target,
        lon_target,
        lat_target,
        target_valid,
    )


def test_completes_missing_drift_over_operational_ice_within_bound():
    values = completion_case()
    u_out, v_out, stats = MODULE.complete_drift_where_ice(
        *values,
        concentration_threshold_percent=15.0,
        maximum_distance_m=30_000.0,
        donor_domain="non_land",
    )

    assert np.isclose(u_out[1, 1], 2.0)
    assert np.isclose(v_out[1, 1], -1.0)
    assert stats["missing_cells"] == 1
    assert stats["repaired_cells"] == 1
    assert stats["unresolved_cells"] == 0
    np.testing.assert_allclose(stats["repair_distances_m"], [25_000.0])


def test_reports_ice_drift_gap_beyond_bound_as_unresolved():
    values = completion_case()
    u_out, v_out, stats = MODULE.complete_drift_where_ice(
        *values,
        concentration_threshold_percent=15.0,
        maximum_distance_m=10_000.0,
        donor_domain="non_land",
    )

    assert np.isnan(u_out[1, 1])
    assert np.isnan(v_out[1, 1])
    assert stats["repaired_cells"] == 0
    assert stats["unresolved_cells"] == 1
    assert np.isclose(stats["minimum_unresolved_distance_m"], 25_000.0)


def test_does_not_require_drift_on_land():
    values = list(completion_case())
    values[-1][1, 1] = False
    _, _, stats = MODULE.complete_drift_where_ice(
        *values,
        concentration_threshold_percent=15.0,
        maximum_distance_m=30_000.0,
        donor_domain="non_land",
    )

    assert stats["ice_cells"] == 0
    assert stats["missing_cells"] == 0
    assert stats["unresolved_cells"] == 0


def test_exclusion_preserves_model_concentration_and_masks_age_field():
    concentration = np.array([[80.0, 5.0], [20.0, 90.0]])
    u_target = np.array([[1.0, np.nan], [np.nan, np.nan]])
    v_target = np.array([[2.0, np.nan], [np.nan, np.nan]])
    target_valid = np.array([[True, True], [True, False]])

    age_sic, supported, excluded = MODULE.make_age_domain_fields(
        concentration,
        u_target,
        v_target,
        target_valid,
        concentration_threshold_percent=10.0,
        unresolved_policy="exclude",
    )

    np.testing.assert_array_equal(
        supported,
        [[True, False], [False, False]],
    )
    np.testing.assert_array_equal(
        excluded,
        [[False, False], [True, False]],
    )
    assert concentration[1, 0] == 20.0
    assert np.isnan(age_sic[1, 0])
    assert age_sic[0, 1] == 5.0
