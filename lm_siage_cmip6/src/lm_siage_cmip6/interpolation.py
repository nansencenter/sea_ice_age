from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def sample_sic_to_elements(
    x_grid_km: np.ndarray,
    y_grid_km: np.ndarray,
    concentration_percent: np.ndarray,
    land_mask: np.ndarray,
    points_yx_km: np.ndarray,
    mode: str = "legacy_zero_fill",
) -> tuple[np.ndarray, np.ndarray]:
    """Sample remapped SIC at element centroids.

    ``model_linear`` uses only finite ocean grid cells and renormalizes the
    bilinear weights. This removes the observational coastal-gap treatment and
    avoids treating land/fill values as zero-concentration ocean.
    ``legacy_zero_fill`` reproduces the previous model-script behavior.
    """
    if mode not in {"model_linear", "legacy_zero_fill"}:
        raise ValueError(f"Unsupported SIC interpolation mode: {mode}")

    concentration = np.ma.asarray(concentration_percent, dtype=np.float64)
    values = np.ma.filled(concentration, np.nan)
    land = np.asarray(land_mask, dtype=bool)
    valid = (
        ~np.ma.getmaskarray(concentration)
        & np.isfinite(values)
        & ~land
    )
    filled = np.where(valid, values, 0.0)
    numerator = RegularGridInterpolator(
        (y_grid_km, x_grid_km), filled,
        method="linear", bounds_error=False, fill_value=0.0,
    )(points_yx_km)
    coverage = RegularGridInterpolator(
        (y_grid_km, x_grid_km), valid.astype(float),
        method="linear", bounds_error=False, fill_value=0.0,
    )(points_yx_km)
    land_at_points = RegularGridInterpolator(
        (y_grid_km, x_grid_km), land.astype(float),
        method="nearest", bounds_error=False, fill_value=1.0,
    )(points_yx_km) > 0.5

    if mode == "model_linear":
        good = coverage > np.finfo(np.float64).eps
        sic = np.full(len(points_yx_km), np.nan, dtype=np.float64)
        sic[good] = numerator[good] / coverage[good]
        missing = ~good & ~land_at_points
    else:
        sic = numerator
        missing = (coverage < 0.5) & ~land_at_points
    sic[land_at_points] = 0.0
    return sic, missing
