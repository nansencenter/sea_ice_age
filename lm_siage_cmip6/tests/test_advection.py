import numpy as np

from lm_siage_cmip6.advection import (
    advect_nodes_carefully,
    boundary_nodes,
    compute_element_mapping,
    fill_missing_velocity_nearest,
    orient_counterclockwise,
    sample_mapping_points,
    sample_velocity_nodes,
    signed_area,
)


def square_mesh():
    x = np.array([0.0, 1.0, 0.0, 1.0])
    y = np.array([0.0, 0.0, 1.0, 1.0])
    triangles = np.array([[0, 1, 2], [1, 3, 2]], dtype=np.int32)
    return x, y, triangles


def test_uniform_translation_preserves_area():
    x, y, triangles = square_mesh()
    triangles = orient_counterclockwise(x, y, triangles)
    x_new, y_new, iterations, affected = advect_nodes_carefully(
        x, y, np.ones(4) * 0.2, np.ones(4) * -0.1, triangles
    )
    np.testing.assert_allclose(x_new, x + 0.2)
    np.testing.assert_allclose(y_new, y - 0.1)
    np.testing.assert_allclose(
        signed_area(x_new, y_new, triangles), signed_area(x, y, triangles)
    )
    assert iterations == 0
    assert affected == 0


def test_identical_mesh_mapping_weights_sum_to_one():
    x, y, triangles = square_mesh()
    grid = np.array([0.25, 0.75])
    domain = np.ones((2, 2), dtype=bool)
    sample_x, sample_y = sample_mapping_points(grid, grid, domain, 3)
    mapping, weights = compute_element_mapping(
        x, y, triangles, x, y, triangles, sample_x, sample_y
    )
    totals = np.zeros(len(triangles))
    np.add.at(totals, mapping[:, 1], weights)
    np.testing.assert_allclose(totals, 1.0)


def test_boundary_nodes_on_square_mesh():
    _, _, triangles = square_mesh()
    np.testing.assert_array_equal(boundary_nodes(triangles), [0, 1, 2, 3])


def test_model_drift_interpolation_does_not_mix_valid_drift_with_missing_zero():
    x_grid = np.array([0.0, 1.0])
    y_grid = np.array([0.0, 1.0])
    u = np.ma.array(
        [[10.0, 10.0], [10.0, 10.0]],
        mask=[[False, True], [False, True]],
    )
    v = np.ma.array(
        [[2.0, 2.0], [2.0, 2.0]],
        mask=[[False, True], [False, True]],
    )
    sampled_u, sampled_v, missing = sample_velocity_nodes(
        x_grid, y_grid, u, v,
        np.array([0.25]), np.array([0.5]),
        interpolation_mode="model_linear", missing_policy="error",
    )
    np.testing.assert_allclose(sampled_u, [10.0])
    np.testing.assert_allclose(sampled_v, [2.0])
    assert not missing.any()


def test_strict_model_drift_policy_leaves_uncovered_nodes_missing():
    values = np.ma.masked_all((2, 2))
    sampled_u, sampled_v, missing = sample_velocity_nodes(
        np.array([0.0, 1.0]), np.array([0.0, 1.0]), values, values,
        np.array([0.5]), np.array([0.5]),
        interpolation_mode="model_linear", missing_policy="error",
    )
    assert missing.all()
    assert np.isnan(sampled_u).all()
    assert np.isnan(sampled_v).all()


def test_nearest_drift_repair_is_bounded():
    u_grid = np.ma.array([[1.0, 2.0]], mask=[[False, False]])
    v_grid = np.ma.array([[3.0, 4.0]], mask=[[False, False]])
    u, v, unresolved, distances = fill_missing_velocity_nearest(
        np.array([0.0, 25.0]), np.array([0.0]), u_grid, v_grid,
        np.array([24.0, 100.0]), np.array([0.0, 0.0]),
        np.array([np.nan, np.nan]), np.array([np.nan, np.nan]),
        np.array([True, True]), maximum_distance_km=50.0,
    )
    np.testing.assert_allclose(u[0], 2.0)
    np.testing.assert_allclose(v[0], 4.0)
    assert not unresolved[0]
    assert unresolved[1]
    np.testing.assert_allclose(distances, [1.0])
