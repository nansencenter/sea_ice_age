from __future__ import annotations

import numpy as np
from pyproj import CRS, Geod, Transformer


def east_north_to_projected(
    longitude: np.ndarray,
    latitude: np.ndarray,
    eastward: np.ndarray,
    northward: np.ndarray,
    target_crs: str | CRS,
    basis_distance_m: float = 10.0,
) -> tuple[np.ma.MaskedArray, np.ma.MaskedArray]:
    """Convert Earth-relative east/north velocity to projected x/y velocity.

    Local projected basis vectors are evaluated geodesically. The result is in
    projected metres per second when the inputs are in metres per second.
    """
    inputs = [np.ma.asarray(value, dtype=np.float64) for value in (
        longitude, latitude, eastward, northward
    )]
    lon, lat, u, v = np.broadcast_arrays(
        *[np.ma.filled(value, np.nan) for value in inputs]
    )
    masks = np.broadcast_arrays(*[np.ma.getmaskarray(value) for value in inputs])
    mask = np.logical_or.reduce(masks) | ~np.isfinite(lon) | ~np.isfinite(lat) | ~np.isfinite(u) | ~np.isfinite(v)
    transformer = Transformer.from_crs("EPSG:4326", target_crs, always_xy=True)
    geod = Geod(ellps="WGS84")
    shape = lon.shape

    def flattened(values: np.ndarray) -> list[float]:
        return np.asarray(values, dtype=np.float64).ravel().tolist()

    def shaped(values: list[float] | np.ndarray) -> np.ndarray:
        return np.asarray(values, dtype=np.float64).reshape(shape)

    lon_flat = flattened(lon)
    lat_flat = flattened(lat)
    x0, y0 = map(shaped, transformer.transform(lon_flat, lat_flat))
    east_azimuth = np.full(lon.size, 90.0).tolist()
    north_azimuth = np.zeros(lon.size).tolist()
    distances = np.full(lon.size, basis_distance_m).tolist()
    lon_e, lat_e, _ = geod.fwd(lon_flat, lat_flat, east_azimuth, distances)
    lon_n, lat_n, _ = geod.fwd(lon_flat, lat_flat, north_azimuth, distances)
    xe, ye = map(shaped, transformer.transform(lon_e, lat_e))
    xn, yn = map(shaped, transformer.transform(lon_n, lat_n))

    bx_e = (xe - x0) / basis_distance_m
    by_e = (ye - y0) / basis_distance_m
    bx_n = (xn - x0) / basis_distance_m
    by_n = (yn - y0) / basis_distance_m
    dx_dt = np.asarray(u) * bx_e + np.asarray(v) * bx_n
    dy_dt = np.asarray(u) * by_e + np.asarray(v) * by_n
    return np.ma.array(dx_dt, mask=mask), np.ma.array(dy_dt, mask=mask)
