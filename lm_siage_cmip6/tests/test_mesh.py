import numpy as np

from lm_siage_cmip6.mesh import build_domain_mask, mesh_quality, triangulate_cell_centres


def test_domain_excludes_land_and_low_latitudes():
    concentration = np.zeros((5, 5))
    concentration[2, 2] = 1.0
    land = np.zeros((5, 5), dtype=bool)
    land[2, 3] = True
    latitude = np.full((5, 5), 80.0)
    latitude[0] = 40.0
    domain = build_domain_mask(concentration, land, latitude, 0.15, 1, 45.0)
    assert domain[2, 2]
    assert not domain[2, 3]
    assert not domain[0].any()


def test_regular_grid_triangulation_has_valid_quality():
    coordinates = np.arange(4, dtype=float) * 25.0
    domain = np.ones((4, 4), dtype=bool)
    x, y, triangles = triangulate_cell_centres(coordinates, coordinates, domain)
    quality = mesh_quality(x, y, triangles)
    assert quality.vertices == 16
    assert quality.triangles == 18
    assert quality.minimum_area_km2 > 0
    assert quality.minimum_quality > 0
