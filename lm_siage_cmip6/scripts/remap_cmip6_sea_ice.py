#!/usr/bin/env python3
"""
remap_cmip6_sea_ice.py

Generic CMIP6 sea-ice remapping workflow.

Inputs:
    siu(time,i,j)      CMIP6 native-grid sea-ice x-velocity
    siv(time,i,j)      CMIP6 native-grid sea-ice y-velocity
    siconc(time,i,j)   CMIP6 native-grid sea-ice concentration
    Grid_example.nc    target LAEA grid

Outputs:
    siu(time,y,x)      eastward sea-ice velocity on LAEA grid
    siv(time,y,x)      northward sea-ice velocity on LAEA grid
    siconc(time,y,x)   sea-ice concentration on LAEA grid

Method:
    1. Estimate native-grid orientation using the perturbation/Jacobian method.
    2. Rotate native siu/siv to geographic east/north on the source grid.
    3. Regrid drift as a true vector:
          east/north -> Earth-centred Cartesian tangent vector
          interpolate Cartesian vector components
          project back to local target tangent plane
          Cartesian tangent vector -> east/north
    4. Regrid siconc as a scalar.
    5. Optionally repair the known pole/origin index block by discarding it and
       filling siu/siv using Cartesian tangent-vector interpolation from the
       surrounding target-grid cells.
    6. Optionally complete missing target-grid drift over operational sea ice
       from a bounded nearby valid target-grid donor, preserving vector
       orientation in geographic coordinates.

No xarray is used.

Example:

python remap_cmip6_sea_ice.py \
  --siu siu_SIday_MODEL_historical_VARIANT_GRID_DATES.nc \
  --siv siv_SIday_MODEL_historical_VARIANT_GRID_DATES.nc \
  --siconc siconc_SIday_MODEL_historical_VARIANT_GRID_DATES.nc \
  --grid ../Grid_example.nc \
  --out seaice_MODEL_YEAR_LAEA.nc \
  --weights weights_MODEL_GRID.npz \
  --time_chunk 50 \
  --write_diagnostics rotation_basis_operational.nc \
  --fix_index_block \
  --index_block_start 214 \
  --index_block_end 218 \
  --index_block_pad 3 \
  --block_method idw \
  --block_idw_power 2.0 \
  --auto_clip_siconc
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime, timezone
from typing import Tuple

import numpy as np
from netCDF4 import Dataset
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
from pyproj import CRS, Transformer, Geod


FILL_VALUE = np.float32(1.0e20)
TARGET_GRID_MAPPING = "Lambert_Azimuthal_Equal_Area"


# =============================================================================
# NetCDF helpers
# =============================================================================

def nc_to_array(var) -> np.ndarray:
    """Read a NetCDF variable as float64, replacing fill/missing values with NaN."""
    arr = var[:]
    if hasattr(arr, "mask"):
        arr = np.asarray(arr.filled(np.nan), dtype=np.float64)
    else:
        arr = np.asarray(arr, dtype=np.float64)

    for att in ["_FillValue", "missing_value"]:
        if att in var.ncattrs():
            try:
                fv = float(np.asarray(var.getncattr(att)).ravel()[0])
                atol = max(abs(fv) * 1.0e-6, 1.0e-30)
                arr[np.isclose(arr, fv, rtol=0.0, atol=atol)] = np.nan
            except Exception:
                pass

    arr[~np.isfinite(arr)] = np.nan
    return arr


def read_time_block_as_float(var, t0: int, count: int) -> np.ndarray:
    """Read var[t0:t0+count,:,:] as float64 with fill/missing converted to NaN."""
    arr = var[t0:t0 + count, :, :]
    if hasattr(arr, "mask"):
        arr = np.asarray(arr.filled(np.nan), dtype=np.float64)
    else:
        arr = np.asarray(arr, dtype=np.float64)

    for att in ["_FillValue", "missing_value"]:
        if att in var.ncattrs():
            try:
                fv = float(np.asarray(var.getncattr(att)).ravel()[0])
                atol = max(abs(fv) * 1.0e-6, 1.0e-30)
                arr[np.isclose(arr, fv, rtol=0.0, atol=atol)] = np.nan
            except Exception:
                pass

    arr[~np.isfinite(arr)] = np.nan
    return arr


def copy_var_attrs(src_var, dst_var, skip=None) -> None:
    if skip is None:
        skip = set()
    for att in src_var.ncattrs():
        if att in skip:
            continue
        try:
            dst_var.setncattr(att, src_var.getncattr(att))
        except Exception:
            pass


def copy_global_attrs(src_ds, dst_ds, skip=None) -> None:
    if skip is None:
        skip = set()
    for att in src_ds.ncattrs():
        if att in skip:
            continue
        try:
            dst_ds.setncattr(att, src_ds.getncattr(att))
        except Exception:
            pass


def get_units(var, fallback: str) -> str:
    return var.getncattr("units") if "units" in var.ncattrs() else fallback


# =============================================================================
# Perturbation/Jacobian native-grid rotation
# =============================================================================

def weighted_lstsq(X: np.ndarray, y: np.ndarray, w: np.ndarray):
    good = np.isfinite(X).all(axis=1) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    if np.count_nonzero(good) < 2:
        return None
    Xg = X[good, :]
    yg = y[good]
    wg = np.sqrt(w[good])
    Xw = Xg * wg[:, None]
    yw = yg * wg
    try:
        beta, _, _, _ = np.linalg.lstsq(Xw, yw, rcond=None)
        return beta
    except Exception:
        return None


def robust_local_fit(
    di: np.ndarray,
    dj: np.ndarray,
    east: np.ndarray,
    north: np.ndarray,
    dist: np.ndarray,
    min_points: int,
    huber: bool = True,
) -> Tuple[float, float, float, float, int, float]:
    valid = (
        np.isfinite(di) & np.isfinite(dj) &
        np.isfinite(east) & np.isfinite(north) &
        np.isfinite(dist) & (dist > 1.0)
    )

    if np.count_nonzero(valid) < min_points:
        return np.nan, np.nan, np.nan, np.nan, int(np.count_nonzero(valid)), np.nan

    di = di[valid]
    dj = dj[valid]
    east = east[valid]
    north = north[valid]

    X = np.column_stack((di, dj))
    index_dist = np.sqrt(di * di + dj * dj)
    w = 1.0 / np.maximum(index_dist, 1.0) ** 2

    beta_e = weighted_lstsq(X, east, w)
    beta_n = weighted_lstsq(X, north, w)

    if beta_e is None or beta_n is None:
        return np.nan, np.nan, np.nan, np.nan, int(len(di)), np.nan

    if huber and len(di) >= max(min_points + 2, 6):
        pred_e = X @ beta_e
        pred_n = X @ beta_n
        res = np.hypot(east - pred_e, north - pred_n)
        med = np.nanmedian(res)
        mad = np.nanmedian(np.abs(res - med))
        scale = 1.4826 * mad if mad > 0 else np.nanmedian(res)
        if np.isfinite(scale) and scale > 0:
            keep = res <= (med + 3.0 * scale)
            if np.count_nonzero(keep) >= min_points:
                beta_e2 = weighted_lstsq(X[keep, :], east[keep], w[keep])
                beta_n2 = weighted_lstsq(X[keep, :], north[keep], w[keep])
                if beta_e2 is not None and beta_n2 is not None:
                    beta_e, beta_n = beta_e2, beta_n2
                    X = X[keep, :]
                    east = east[keep]
                    north = north[keep]

    # Columns of the local Jacobian:
    # +i vector = (d east / d i, d north / d i)
    # +j vector = (d east / d j, d north / d j)
    vi_e = beta_e[0]
    vj_e = beta_e[1]
    vi_n = beta_n[0]
    vj_n = beta_n[1]

    ni = np.hypot(vi_e, vi_n)
    nj = np.hypot(vj_e, vj_n)

    if not np.isfinite(ni) or ni <= 0 or not np.isfinite(nj) or nj <= 0:
        return np.nan, np.nan, np.nan, np.nan, int(len(di)), np.nan

    bi_e = vi_e / ni
    bi_n = vi_n / ni
    bj_e = vj_e / nj
    bj_n = vj_n / nj

    pred_e = X @ beta_e
    pred_n = X @ beta_n
    rms = np.sqrt(np.nanmean((east - pred_e) ** 2 + (north - pred_n) ** 2))

    return bi_e, bi_n, bj_e, bj_n, int(len(di)), rms


def compute_perturbation_basis(
    lon: np.ndarray,
    lat: np.ndarray,
    stencil_radius: int = 2,
    max_neighbor_dist: float = 400000.0,
    min_points: int = 5,
    periodic_j: bool = True,
    robust: bool = True,
) -> Tuple[np.ndarray, ...]:
    lon = np.asarray(lon, dtype=np.float64)
    lat = np.asarray(lat, dtype=np.float64)
    ni, nj = lon.shape
    geod = Geod(ellps="WGS84")

    bi_e = np.full((ni, nj), np.nan, dtype=np.float64)
    bi_n = np.full((ni, nj), np.nan, dtype=np.float64)
    bj_e = np.full((ni, nj), np.nan, dtype=np.float64)
    bj_n = np.full((ni, nj), np.nan, dtype=np.float64)
    n_used = np.zeros((ni, nj), dtype=np.int16)
    rms = np.full((ni, nj), np.nan, dtype=np.float64)

    offsets = []
    for ddi in range(-stencil_radius, stencil_radius + 1):
        for ddj in range(-stencil_radius, stencil_radius + 1):
            if ddi == 0 and ddj == 0:
                continue
            if np.sqrt(ddi * ddi + ddj * ddj) <= stencil_radius + 1.0e-12:
                offsets.append((ddi, ddj))

    for ii in range(ni):
        if ii % 25 == 0:
            print(f"  basis row {ii}/{ni - 1}")

        for jj in range(nj):
            lon0 = lon[ii, jj]
            lat0 = lat[ii, jj]

            if not (np.isfinite(lon0) and np.isfinite(lat0)):
                continue

            di_list = []
            dj_list = []
            east_list = []
            north_list = []
            dist_list = []

            for ddi, ddj in offsets:
                ii2 = ii + ddi
                jj2 = jj + ddj

                if ii2 < 0 or ii2 >= ni:
                    continue

                if periodic_j:
                    jj2 = jj2 % nj
                    dj_fit = ddj
                else:
                    if jj2 < 0 or jj2 >= nj:
                        continue
                    dj_fit = ddj

                lon1 = lon[ii2, jj2]
                lat1 = lat[ii2, jj2]

                if not (np.isfinite(lon1) and np.isfinite(lat1)):
                    continue

                try:
                    az, _, dist = geod.inv(lon0, lat0, lon1, lat1)
                except Exception:
                    continue

                if not np.isfinite(dist) or dist <= 1.0 or dist > max_neighbor_dist:
                    continue

                a = np.deg2rad(az)
                east = dist * np.sin(a)
                north = dist * np.cos(a)

                di_list.append(float(ddi))
                dj_list.append(float(dj_fit))
                east_list.append(east)
                north_list.append(north)
                dist_list.append(dist)

            out = robust_local_fit(
                np.asarray(di_list),
                np.asarray(dj_list),
                np.asarray(east_list),
                np.asarray(north_list),
                np.asarray(dist_list),
                min_points=min_points,
                huber=robust,
            )

            bi_e[ii, jj], bi_n[ii, jj], bj_e[ii, jj], bj_n[ii, jj], n_used[ii, jj], rms[ii, jj] = out

    return bi_e, bi_n, bj_e, bj_n, n_used, rms


def fill_basis_nearest(
    bi_e: np.ndarray,
    bi_n: np.ndarray,
    bj_e: np.ndarray,
    bj_n: np.ndarray,
    max_iter: int = 10,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    bi_e = bi_e.copy()
    bi_n = bi_n.copy()
    bj_e = bj_e.copy()
    bj_n = bj_n.copy()

    filled_total = 0
    ni, nj = bi_e.shape

    for _ in range(max_iter):
        missing = ~(np.isfinite(bi_e) & np.isfinite(bi_n) & np.isfinite(bj_e) & np.isfinite(bj_n))
        if not np.any(missing):
            break

        to_fill = np.argwhere(missing)
        filled_this = 0

        bi_e_new = bi_e.copy()
        bi_n_new = bi_n.copy()
        bj_e_new = bj_e.copy()
        bj_n_new = bj_n.copy()

        for ii, jj in to_fill:
            vals_bi_e = []
            vals_bi_n = []
            vals_bj_e = []
            vals_bj_n = []

            for ddi in [-1, 0, 1]:
                for ddj in [-1, 0, 1]:
                    if ddi == 0 and ddj == 0:
                        continue
                    ii2 = ii + ddi
                    jj2 = (jj + ddj) % nj
                    if ii2 < 0 or ii2 >= ni:
                        continue
                    good = (
                        np.isfinite(bi_e[ii2, jj2]) and np.isfinite(bi_n[ii2, jj2]) and
                        np.isfinite(bj_e[ii2, jj2]) and np.isfinite(bj_n[ii2, jj2])
                    )
                    if good:
                        vals_bi_e.append(bi_e[ii2, jj2])
                        vals_bi_n.append(bi_n[ii2, jj2])
                        vals_bj_e.append(bj_e[ii2, jj2])
                        vals_bj_n.append(bj_n[ii2, jj2])

            if len(vals_bi_e) >= 1:
                vie = np.nanmean(vals_bi_e)
                vin = np.nanmean(vals_bi_n)
                vje = np.nanmean(vals_bj_e)
                vjn = np.nanmean(vals_bj_n)

                ni_norm = np.hypot(vie, vin)
                nj_norm = np.hypot(vje, vjn)

                if ni_norm > 0 and nj_norm > 0:
                    bi_e_new[ii, jj] = vie / ni_norm
                    bi_n_new[ii, jj] = vin / ni_norm
                    bj_e_new[ii, jj] = vje / nj_norm
                    bj_n_new[ii, jj] = vjn / nj_norm
                    filled_this += 1

        bi_e, bi_n, bj_e, bj_n = bi_e_new, bi_n_new, bj_e_new, bj_n_new
        filled_total += filled_this

        if filled_this == 0:
            break

    return bi_e, bi_n, bj_e, bj_n, filled_total


def rotate_native_to_east_north(
    siu: np.ndarray,
    siv: np.ndarray,
    bi_e: np.ndarray,
    bi_n: np.ndarray,
    bj_e: np.ndarray,
    bj_n: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    # siu positive along native +j; siv positive along native +i.
    u_east = siu * bj_e + siv * bi_e
    v_north = siu * bj_n + siv * bi_n
    return u_east, v_north


# =============================================================================
# Vector basis conversion: east/north <-> 3-D Cartesian tangent vectors
# =============================================================================

def local_basis_vectors(lon_deg, lat_deg):
    lon = np.deg2rad(lon_deg)
    lat = np.deg2rad(lat_deg)

    sinlon = np.sin(lon)
    coslon = np.cos(lon)
    sinlat = np.sin(lat)
    coslat = np.cos(lat)

    east = np.stack((-sinlon, coslon, np.zeros_like(lon)), axis=-1)
    north = np.stack((-sinlat * coslon, -sinlat * sinlon, coslat), axis=-1)
    radial = np.stack((coslat * coslon, coslat * sinlon, sinlat), axis=-1)

    return east, north, radial


def en_to_cartesian(u_east, v_north, lon_deg, lat_deg):
    east, north, _ = local_basis_vectors(lon_deg, lat_deg)
    return u_east[..., None] * east + v_north[..., None] * north


def cartesian_to_en(vec_cart, lon_deg, lat_deg):
    east, north, radial = local_basis_vectors(lon_deg, lat_deg)

    radial_component = np.sum(vec_cart * radial, axis=-1)
    vec_tangent = vec_cart - radial_component[..., None] * radial

    u = np.sum(vec_tangent * east, axis=-1)
    v = np.sum(vec_tangent * north, axis=-1)
    return u, v


# =============================================================================
# Regridding helpers
# =============================================================================

def make_target_crs(grid_ds) -> CRS:
    if TARGET_GRID_MAPPING not in grid_ds.variables:
        raise KeyError(f"Target grid is missing {TARGET_GRID_MAPPING}")
    gm = grid_ds.variables[TARGET_GRID_MAPPING]
    if "proj4_string" in gm.ncattrs():
        return CRS.from_proj4(gm.getncattr("proj4_string"))
    return CRS.from_cf({att: gm.getncattr(att) for att in gm.ncattrs()})


def get_target_valid_mask(
    grid_ds,
    lon_tgt: np.ndarray,
    lat_tgt: np.ndarray,
    mode: str = "land_only",
    mask_var: str = "status_flag",
    land_value: int = 1,
    nominal_value: int = 0,
) -> np.ndarray:
    finite = np.isfinite(lon_tgt) & np.isfinite(lat_tgt)

    if mode == "none":
        return finite

    if mask_var not in grid_ds.variables:
        print(f"Warning: target mask variable '{mask_var}' not found. Using finite lon/lat only.")
        return finite

    flag = nc_to_array(grid_ds.variables[mask_var])
    if flag.ndim == 3:
        flag = flag[0, :, :]
    elif flag.ndim != 2:
        raise ValueError(f"Unsupported target mask dimensions for {mask_var}: {flag.shape}")

    if flag.shape != lon_tgt.shape:
        raise ValueError(f"Target mask shape {flag.shape} does not match target lon/lat {lon_tgt.shape}")

    if mode == "land_only":
        return finite & np.isfinite(flag) & (flag != land_value)

    if mode == "nominal":
        return finite & np.isfinite(flag) & (flag == nominal_value)

    raise ValueError(f"Unsupported target_mask_mode: {mode}")


def build_weights(
    x_src_flat: np.ndarray,
    y_src_flat: np.ndarray,
    source_good: np.ndarray,
    x_tgt_flat: np.ndarray,
    y_tgt_flat: np.ndarray,
    target_good: np.ndarray,
    k: int,
    max_dist: float,
    power: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    source_query_good = source_good & np.isfinite(x_src_flat) & np.isfinite(y_src_flat)
    target_query_good = target_good & np.isfinite(x_tgt_flat) & np.isfinite(y_tgt_flat)

    n_source = int(np.count_nonzero(source_query_good))
    n_target = int(np.count_nonzero(target_query_good))

    if n_source < k:
        raise ValueError(f"Not enough valid source points for k={k}: {n_source}")
    if n_target == 0:
        raise ValueError("No valid finite target points available for interpolation")

    print(f"Valid source points for weights: {n_source}")
    print(f"Valid target points for weights: {n_target}")

    source_indices = np.where(source_query_good)[0]
    target_indices = np.where(target_query_good)[0]

    tree = cKDTree(np.column_stack((x_src_flat[source_query_good], y_src_flat[source_query_good])))
    target_points = np.column_stack((x_tgt_flat[target_query_good], y_tgt_flat[target_query_good]))

    dist, idx_local = tree.query(target_points, k=k)

    if k == 1:
        dist = dist[:, None]
        idx_local = idx_local[:, None]

    idx_global = source_indices[idx_local]

    w = 1.0 / np.maximum(dist, 1.0) ** power
    w = w / np.sum(w, axis=1)[:, None]

    too_far = dist[:, 0] > max_dist
    w[too_far, :] = np.nan

    full_indices = np.zeros((x_tgt_flat.size, k), dtype=np.int64)
    full_weights = np.full((x_tgt_flat.size, k), np.nan, dtype=np.float64)
    full_nearest = np.full(x_tgt_flat.size, -1, dtype=np.int64)

    full_indices[target_indices, :] = idx_global
    full_weights[target_indices, :] = w
    full_nearest[target_indices] = idx_global[:, 0]
    full_nearest[target_indices[too_far]] = -1

    return full_indices, full_weights, full_nearest


def apply_weights(source_flat: np.ndarray, indices: np.ndarray, weights: np.ndarray) -> np.ndarray:
    source_flat = np.asarray(source_flat, dtype=np.float64)
    values = source_flat[indices]

    valid = np.isfinite(values) & np.isfinite(weights)
    numerator = np.nansum(np.where(valid, values * weights, np.nan), axis=1)
    denominator = np.nansum(np.where(valid, weights, np.nan), axis=1)

    with np.errstate(invalid="ignore", divide="ignore"):
        out = numerator / denominator

    out[denominator <= 0.0] = np.nan
    out[~np.isfinite(weights).any(axis=1)] = np.nan
    return out


def vector_cartesian_remap(
    u_src: np.ndarray,
    v_src: np.ndarray,
    lon_src: np.ndarray,
    lat_src: np.ndarray,
    lon_tgt: np.ndarray,
    lat_tgt: np.ndarray,
    indices: np.ndarray,
    weights: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Remap an east/north vector field by converting source vectors to 3-D
    Cartesian tangent vectors, applying interpolation weights to the Cartesian
    components, and projecting back to local east/north on the target grid.
    """
    valid_src = (
        np.isfinite(u_src) & np.isfinite(v_src) &
        np.isfinite(lon_src) & np.isfinite(lat_src)
    )

    vec_src = np.full(lon_src.shape + (3,), np.nan, dtype=np.float64)
    vec_src[valid_src, :] = en_to_cartesian(
        u_src[valid_src],
        v_src[valid_src],
        lon_src[valid_src],
        lat_src[valid_src],
    )

    comp_tgt = []
    for cc in range(3):
        comp = apply_weights(vec_src[..., cc].ravel(), indices, weights)
        comp_tgt.append(comp.reshape(lon_tgt.shape))

    vec_tgt = np.stack(comp_tgt, axis=-1)
    u_tgt, v_tgt = cartesian_to_en(vec_tgt, lon_tgt, lat_tgt)

    bad = ~np.isfinite(vec_tgt).all(axis=-1)
    u_tgt[bad] = np.nan
    v_tgt[bad] = np.nan

    return u_tgt, v_tgt


def make_drift_extent_mask(
    u_src_flat: np.ndarray,
    v_src_flat: np.ndarray,
    indices: np.ndarray,
    weights: np.ndarray,
    nearest_index: np.ndarray,
    mode: str,
    min_valid_weight: float,
) -> np.ndarray:
    finite_source = np.isfinite(u_src_flat) & np.isfinite(v_src_flat)

    if mode == "none":
        return np.ones(nearest_index.shape, dtype=bool)

    if mode == "nearest":
        out = np.zeros(nearest_index.shape, dtype=bool)
        ok = nearest_index >= 0
        out[ok] = finite_source[nearest_index[ok]]
        return out

    if mode == "weighted_fraction":
        valid_neighbour = finite_source[indices] & np.isfinite(weights)
        valid_weight = np.nansum(np.where(valid_neighbour, weights, 0.0), axis=1)
        return valid_weight >= min_valid_weight

    if mode == "all":
        valid_neighbour = finite_source[indices] & np.isfinite(weights)
        return np.all(valid_neighbour, axis=1)

    if mode == "any":
        valid_neighbour = finite_source[indices] & np.isfinite(weights)
        return np.any(valid_neighbour, axis=1)

    raise ValueError(f"Unsupported drift_extent_mode: {mode}")


def complete_drift_where_ice(
    u_target: np.ndarray,
    v_target: np.ndarray,
    concentration_target: np.ndarray,
    x_target: np.ndarray,
    y_target: np.ndarray,
    lon_target: np.ndarray,
    lat_target: np.ndarray,
    target_valid: np.ndarray,
    concentration_threshold_percent: float,
    maximum_distance_m: float,
    donor_domain: str = "non_land",
) -> Tuple[np.ndarray, np.ndarray, dict]:
    """Fill missing drift over operational ice from bounded target-grid donors.

    Donor vectors are converted to Earth-centred Cartesian tangent vectors and
    projected onto the recipient cell's local east/north basis. This avoids
    treating east/north components at separated geographic points as if their
    local coordinate bases were identical.
    """
    if maximum_distance_m <= 0:
        raise ValueError("drift fill maximum distance must be positive")
    if not 0.0 <= concentration_threshold_percent <= 100.0:
        raise ValueError(
            "operational concentration threshold must be between 0 and 100 percent"
        )
    if donor_domain not in {"non_land", "ice"}:
        raise ValueError(f"Unsupported drift fill donor domain: {donor_domain}")

    u_out = np.asarray(u_target, dtype=np.float64).copy()
    v_out = np.asarray(v_target, dtype=np.float64).copy()
    concentration = np.asarray(concentration_target, dtype=np.float64)
    valid_target = np.asarray(target_valid, dtype=bool)

    operational_ice = (
        valid_target
        & np.isfinite(concentration)
        & (concentration >= concentration_threshold_percent)
    )
    coordinate_valid = (
        np.isfinite(x_target)
        & np.isfinite(y_target)
        & np.isfinite(lon_target)
        & np.isfinite(lat_target)
    )
    valid_drift = (
        valid_target
        & coordinate_valid
        & np.isfinite(u_out)
        & np.isfinite(v_out)
    )
    missing_ice_drift = operational_ice & ~valid_drift

    stats = {
        "ice_cells": int(np.count_nonzero(operational_ice)),
        "missing_cells": int(np.count_nonzero(missing_ice_drift)),
        "repaired_cells": 0,
        "unresolved_cells": int(np.count_nonzero(missing_ice_drift)),
        "repair_distances_m": np.empty(0, dtype=np.float64),
        "minimum_unresolved_distance_m": np.nan,
    }
    if not missing_ice_drift.any():
        return u_out, v_out, stats

    donor_mask = valid_drift.copy()
    if donor_domain == "ice":
        donor_mask &= operational_ice
    if not donor_mask.any():
        return u_out, v_out, stats

    donor_rows, donor_cols = np.where(donor_mask)
    missing_rows, missing_cols = np.where(missing_ice_drift)
    donor_points = np.column_stack((
        x_target[donor_rows, donor_cols],
        y_target[donor_rows, donor_cols],
    ))
    missing_points = np.column_stack((
        x_target[missing_rows, missing_cols],
        y_target[missing_rows, missing_cols],
    ))

    tree = cKDTree(donor_points)
    distances_m, donor_local_indices = tree.query(missing_points, k=1)
    repairable = np.isfinite(distances_m) & (distances_m <= maximum_distance_m)

    if repairable.any():
        recipient_rows = missing_rows[repairable]
        recipient_cols = missing_cols[repairable]
        selected_donors = donor_local_indices[repairable]
        source_rows = donor_rows[selected_donors]
        source_cols = donor_cols[selected_donors]

        donor_vectors = en_to_cartesian(
            u_out[source_rows, source_cols],
            v_out[source_rows, source_cols],
            lon_target[source_rows, source_cols],
            lat_target[source_rows, source_cols],
        )
        repaired_u, repaired_v = cartesian_to_en(
            donor_vectors,
            lon_target[recipient_rows, recipient_cols],
            lat_target[recipient_rows, recipient_cols],
        )
        u_out[recipient_rows, recipient_cols] = repaired_u
        v_out[recipient_rows, recipient_cols] = repaired_v

    unresolved_distances = distances_m[~repairable]
    stats.update({
        "repaired_cells": int(np.count_nonzero(repairable)),
        "unresolved_cells": int(np.count_nonzero(~repairable)),
        "repair_distances_m": np.asarray(distances_m[repairable], dtype=np.float64),
        "minimum_unresolved_distance_m": (
            float(np.min(unresolved_distances))
            if unresolved_distances.size else np.nan
        ),
    })
    return u_out, v_out, stats


def make_age_domain_fields(
    concentration_target: np.ndarray,
    u_target: np.ndarray,
    v_target: np.ndarray,
    target_valid: np.ndarray,
    concentration_threshold_percent: float,
    unresolved_policy: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build operational concentration and explicit target drift support."""
    if unresolved_policy not in {"error", "exclude"}:
        raise ValueError(f"Unsupported unresolved drift policy: {unresolved_policy}")

    concentration = np.asarray(concentration_target, dtype=np.float64)
    valid_target = np.asarray(target_valid, dtype=bool)
    drift_supported = (
        valid_target
        & np.isfinite(u_target)
        & np.isfinite(v_target)
    )
    unresolved_ice = (
        valid_target
        & np.isfinite(concentration)
        & (concentration >= concentration_threshold_percent)
        & ~drift_supported
    )
    concentration_for_age = concentration.copy()
    if unresolved_policy == "exclude":
        concentration_for_age[unresolved_ice] = np.nan
    return concentration_for_age, drift_supported, unresolved_ice


# =============================================================================
# Final pole/origin index-block repair
# =============================================================================

def build_index_block_masks(ny, nx, start, end, pad):
    if start < 0 or end < start:
        raise ValueError("Invalid block range")
    if end >= ny or end >= nx:
        raise IndexError(f"Block end={end} outside shape ny={ny}, nx={nx}")

    y0 = int(start)
    y1 = int(end) + 1
    x0 = int(start)
    x1 = int(end) + 1

    p = max(1, int(pad))

    wy0 = max(0, y0 - p)
    wy1 = min(ny, y1 + p)
    wx0 = max(0, x0 - p)
    wx1 = min(nx, x1 + p)

    inner = np.zeros((ny, nx), dtype=bool)
    window = np.zeros((ny, nx), dtype=bool)

    inner[y0:y1, x0:x1] = True
    window[wy0:wy1, wx0:wx1] = True

    ring = window & (~inner)
    return inner, ring, (y0, y1, x0, x1), (wy0, wy1, wx0, wx1)


def cartesian_vector_block_repair(
    u_field,
    v_field,
    x2d,
    y2d,
    lon2d,
    lat2d,
    inner,
    ring,
    method="idw",
    idw_power=2.0,
    nearest_fallback=True,
):
    """
    Discard the selected target-grid block and refill it using Cartesian
    tangent-vector interpolation from the surrounding target-grid ring.
    """
    u = np.asarray(u_field, dtype=np.float64).copy()
    v = np.asarray(v_field, dtype=np.float64).copy()

    # Throw away the central block first.
    u[inner] = np.nan
    v[inner] = np.nan

    source = (
        ring
        & np.isfinite(u)
        & np.isfinite(v)
        & np.isfinite(x2d)
        & np.isfinite(y2d)
        & np.isfinite(lon2d)
        & np.isfinite(lat2d)
    )

    target = (
        inner
        & np.isfinite(x2d)
        & np.isfinite(y2d)
        & np.isfinite(lon2d)
        & np.isfinite(lat2d)
    )

    nsource = int(np.count_nonzero(source))
    ntarget = int(np.count_nonzero(target))

    if nsource < 1 or ntarget == 0:
        return u, v, 0, nsource, ntarget

    src_x = x2d[source].astype(np.float64)
    src_y = y2d[source].astype(np.float64)
    src_lon = lon2d[source].astype(np.float64)
    src_lat = lat2d[source].astype(np.float64)

    tgt_x = x2d[target].astype(np.float64)
    tgt_y = y2d[target].astype(np.float64)
    tgt_lon = lon2d[target].astype(np.float64)
    tgt_lat = lat2d[target].astype(np.float64)

    src_vec = en_to_cartesian(u[source], v[source], src_lon, src_lat)

    if method == "linear" and nsource >= 3:
        src_points = np.column_stack((src_x, src_y))
        tgt_points = np.column_stack((tgt_x, tgt_y))

        comps = []
        for cc in range(3):
            try:
                vals = griddata(src_points, src_vec[:, cc], tgt_points, method="linear")
            except Exception:
                vals = np.full(ntarget, np.nan, dtype=np.float64)

            if nearest_fallback:
                missing = ~np.isfinite(vals)
                if np.any(missing):
                    try:
                        vals_nn = griddata(src_points, src_vec[:, cc], tgt_points[missing], method="nearest")
                        vals[missing] = vals_nn
                    except Exception:
                        pass
            comps.append(vals)

        tgt_vec = np.column_stack(comps)

    elif method == "idw":
        tgt_vec = np.full((ntarget, 3), np.nan, dtype=np.float64)
        for kk in range(ntarget):
            d = np.hypot(src_x - tgt_x[kk], src_y - tgt_y[kk])
            good = np.isfinite(d) & np.isfinite(src_vec).all(axis=1)

            if not np.any(good):
                continue

            dg = d[good]
            vg = src_vec[good, :]

            if np.any(dg < 1.0e-9):
                tgt_vec[kk, :] = vg[np.argmin(dg), :]
                continue

            w = 1.0 / np.maximum(dg, 1.0) ** float(idw_power)
            wsum = np.sum(w)
            if np.isfinite(wsum) and wsum > 0:
                tgt_vec[kk, :] = np.sum(vg * w[:, None], axis=0) / wsum

    else:
        raise ValueError(f"Unsupported block repair method: {method}")

    ok = np.isfinite(tgt_vec).all(axis=1)
    if np.any(ok):
        tgt_u, tgt_v = cartesian_to_en(tgt_vec[ok, :], tgt_lon[ok], tgt_lat[ok])
        yy, xx = np.where(target)
        u[yy[ok], xx[ok]] = tgt_u
        v[yy[ok], xx[ok]] = tgt_v

    return u, v, int(np.count_nonzero(ok)), nsource, ntarget


# =============================================================================
# Diagnostics and output creation
# =============================================================================

def write_rotation_diagnostics(
    path: str,
    lon: np.ndarray,
    lat: np.ndarray,
    bi_e: np.ndarray,
    bi_n: np.ndarray,
    bj_e: np.ndarray,
    bj_n: np.ndarray,
    n_used: np.ndarray,
    rms: np.ndarray,
) -> None:
    ni, nj = lon.shape
    with Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("i", ni)
        ds.createDimension("j", nj)

        lonv = ds.createVariable("longitude", "f8", ("i", "j"))
        latv = ds.createVariable("latitude", "f8", ("i", "j"))
        lonv[:, :] = lon
        latv[:, :] = lat
        lonv.units = "degrees_east"
        latv.units = "degrees_north"

        diag_vars = {
            "basis_i_east": bi_e,
            "basis_i_north": bi_n,
            "basis_j_east": bj_e,
            "basis_j_north": bj_n,
            "basis_fit_rms_m": rms,
        }

        for name, arr in diag_vars.items():
            v = ds.createVariable(name, "f4", ("i", "j"), fill_value=FILL_VALUE, zlib=True)
            v[:, :] = np.where(np.isfinite(arr), arr, FILL_VALUE).astype(np.float32)

        nv = ds.createVariable("basis_fit_n_points", "i2", ("i", "j"), zlib=True)
        nv[:, :] = n_used.astype(np.int16)

        ds.title = "Rotation basis diagnostics for perturbation/Jacobian method"


def create_output(
    out_path: str,
    siu_path: str,
    siv_path: str,
    siconc_path: str,
    grid_path: str,
    nt: int,
    ny: int,
    nx: int,
    siconc_var_name: str,
    target_mask_mode: str,
) -> Dataset:
    siu_ds = Dataset(siu_path, "r")
    siv_ds = Dataset(siv_path, "r")
    sic_ds = Dataset(siconc_path, "r")
    grid_ds = Dataset(grid_path, "r")
    out = Dataset(out_path, "w", format="NETCDF4")

    try:
        copy_global_attrs(sic_ds, out)

        out.title = "CMIP6 sea-ice drift and concentration remapped to Arctic LAEA grid"
        out.source_siu_file = os.path.basename(siu_path)
        out.source_siv_file = os.path.basename(siv_path)
        out.source_siconc_file = os.path.basename(siconc_path)
        out.target_grid_file = os.path.basename(grid_path)
        out.target_mask_mode = target_mask_mode
        out.rotation_method = "perturbation/Jacobian native-grid basis"
        out.vector_remapping_method = (
            "Source east/north drift was converted to Earth-centred Cartesian tangent vectors; "
            "Cartesian vector components were interpolated; interpolated vectors were projected "
            "onto the target local tangent plane and decomposed back into east/north components."
        )
        out.history = (
            f"{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}: "
            "rotated native siu/siv to east/north, remapped drift with Cartesian tangent-vector "
            "interpolation, remapped siconc as scalar, and wrote siu/siv/siconc to LAEA grid."
        )
        out.grid = "remapped to target Lambert Azimuthal Equal Area grid"
        out.grid_label = "laea"
        out.Conventions = "CF-1.7"

        out.createDimension("time", None)
        if "bnds" in sic_ds.dimensions:
            out.createDimension("bnds", len(sic_ds.dimensions["bnds"]))
        elif "bnds" in siu_ds.dimensions:
            out.createDimension("bnds", len(siu_ds.dimensions["bnds"]))
        else:
            out.createDimension("bnds", 2)

        out.createDimension("y", ny)
        out.createDimension("x", nx)

        time_source = sic_ds if "time" in sic_ds.variables else siu_ds

        if "time" in time_source.variables:
            vin = time_source.variables["time"]
            vout = out.createVariable("time", vin.datatype, ("time",), zlib=True)
            vout[:] = vin[:]
            copy_var_attrs(vin, vout)

        if "time_bnds" in time_source.variables:
            vin = time_source.variables["time_bnds"]
            vout = out.createVariable("time_bnds", vin.datatype, ("time", "bnds"), zlib=True)
            vout[:, :] = vin[:, :]
            copy_var_attrs(vin, vout)

        if "strlen" in sic_ds.dimensions and "type" in sic_ds.variables:
            out.createDimension("strlen", len(sic_ds.dimensions["strlen"]))
            vin = sic_ds.variables["type"]
            vout = out.createVariable("type", vin.datatype, vin.dimensions)
            vout[:] = vin[:]
            copy_var_attrs(vin, vout)

        for coord in ["x", "y", "longitude", "latitude"]:
            vin = grid_ds.variables[coord]
            vout = out.createVariable(coord, vin.datatype, vin.dimensions, zlib=(coord in ["longitude", "latitude"]))
            vout[:] = vin[:]
            copy_var_attrs(vin, vout)

        if "status_flag" in grid_ds.variables:
            vin = grid_ds.variables["status_flag"]
            vout = out.createVariable("status_flag", vin.datatype, vin.dimensions, zlib=True)
            vout[:] = vin[:]
            copy_var_attrs(vin, vout)

        gm_in = grid_ds.variables[TARGET_GRID_MAPPING]
        gm_out = out.createVariable(TARGET_GRID_MAPPING, gm_in.datatype)
        try:
            gm_out[...] = gm_in[...]
        except Exception:
            pass
        copy_var_attrs(gm_in, gm_out)

        chunks = (1, ny, nx)

        siu_out = out.createVariable("siu", "f4", ("time", "y", "x"), zlib=True, complevel=4, fill_value=FILL_VALUE, chunksizes=chunks)
        siv_out = out.createVariable("siv", "f4", ("time", "y", "x"), zlib=True, complevel=4, fill_value=FILL_VALUE, chunksizes=chunks)
        sic_out = out.createVariable(siconc_var_name, "f4", ("time", "y", "x"), zlib=True, complevel=4, fill_value=FILL_VALUE, chunksizes=chunks)
        age_sic_out = out.createVariable(
            "siconc_for_age", "f4", ("time", "y", "x"),
            zlib=True, complevel=4, fill_value=FILL_VALUE, chunksizes=chunks,
        )
        drift_support_out = out.createVariable(
            "drift_supported", "u1", ("time", "y", "x"),
            zlib=True, complevel=4, fill_value=np.uint8(255), chunksizes=chunks,
        )

        siu_out.units = get_units(siu_ds.variables["siu"], "m s-1")
        siv_out.units = get_units(siv_ds.variables["siv"], "m s-1")
        siu_out.standard_name = "sea_ice_eastward_velocity"
        siv_out.standard_name = "sea_ice_northward_velocity"
        siu_out.long_name = "eastward sea-ice velocity on target LAEA grid"
        siv_out.long_name = "northward sea-ice velocity on target LAEA grid"
        siu_out.coordinates = "latitude longitude"
        siv_out.coordinates = "latitude longitude"
        siu_out.grid_mapping = TARGET_GRID_MAPPING
        siv_out.grid_mapping = TARGET_GRID_MAPPING
        siu_out.comment = "Variable name kept as siu by request. Values are eastward velocity after native-grid rotation and Cartesian tangent-vector remapping."
        siv_out.comment = "Variable name kept as siv by request. Values are northward velocity after native-grid rotation and Cartesian tangent-vector remapping."

        sic_in = sic_ds.variables[siconc_var_name]
        copy_var_attrs(
            sic_in,
            sic_out,
            skip={"_FillValue", "missing_value", "coordinates", "grid_mapping", "cell_measures", "_ChunkSizes"},
        )
        sic_out.setncattr("missing_value", FILL_VALUE)
        sic_out.coordinates = "latitude longitude"
        sic_out.grid_mapping = TARGET_GRID_MAPPING
        old_comment = sic_in.getncattr("comment") if "comment" in sic_in.ncattrs() else ""
        sic_out.comment = (
            old_comment
            + " Remapped to target LAEA grid as a scalar; land target cells "
            "kept missing. Original remapped model concentration is preserved "
            "even where drift is unsupported."
        )
        sic_out.ancillary_variables = "drift_supported"
        copy_var_attrs(
            sic_in,
            age_sic_out,
            skip={
                "_FillValue", "missing_value", "coordinates", "grid_mapping",
                "cell_measures", "_ChunkSizes", "comment", "long_name",
            },
        )
        age_sic_out.setncattr("missing_value", FILL_VALUE)
        age_sic_out.coordinates = "latitude longitude"
        age_sic_out.grid_mapping = TARGET_GRID_MAPPING
        age_sic_out.ancillary_variables = "drift_supported"
        age_sic_out.long_name = (
            "sea-ice concentration restricted to the drift-supported age domain"
        )
        age_sic_out.comment = (
            "Copy of remapped siconc used by the sea-ice-age calculation. "
            "Cells at or above the configured drift-completion concentration "
            "threshold are masked when bounded completion cannot provide "
            "finite siu and siv. See siconc for the unmodified remapped model "
            "concentration."
        )
        drift_support_out.long_name = (
            "finite sea-ice drift support on the target grid"
        )
        drift_support_out.flag_values = np.asarray([0, 1], dtype=np.uint8)
        drift_support_out.flag_meanings = "unsupported supported"
        drift_support_out.coordinates = "latitude longitude"
        drift_support_out.grid_mapping = TARGET_GRID_MAPPING
        drift_support_out.comment = (
            "One where finite siu and siv are available after bounded drift "
            "completion; zero where the non-land target grid lacks finite "
            "drift; fill outside the configured target domain."
        )

        siu_ds.close()
        siv_ds.close()
        sic_ds.close()
        grid_ds.close()
        return out

    except Exception:
        out.close()
        siu_ds.close()
        siv_ds.close()
        sic_ds.close()
        grid_ds.close()
        raise


# =============================================================================
# Main
# =============================================================================

def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Remap CMIP6 sea-ice vectors and concentration to LAEA."
    )

    p.add_argument("--siu", required=True, help="Input CMIP6 siu file")
    p.add_argument("--siv", required=True, help="Input CMIP6 siv file")
    p.add_argument("--siconc", required=True, help="Input CMIP6 siconc file")
    p.add_argument("--grid", required=True, help="Target LAEA grid file")
    p.add_argument("--out", required=True, help="Output NetCDF with siu, siv, siconc on LAEA")
    p.add_argument("--siconc_var", default="siconc", help="Name of sea-ice concentration variable")

    # Rotation
    p.add_argument("--stencil_radius", type=int, default=2)
    p.add_argument("--max_neighbor_dist", type=float, default=400000.0)
    p.add_argument("--min_points", type=int, default=5)
    p.add_argument("--no_periodic_j", action="store_true")
    p.add_argument("--no_robust", action="store_true")
    p.add_argument("--no_fill_basis", action="store_true")
    p.add_argument("--write_diagnostics", default=None)

    # Remapping weights
    p.add_argument("--weights", default=None, help="Optional .npz weights file to save/reuse")
    p.add_argument("--k", type=int, default=4, help="Number of nearest neighbours for IDW remapping")
    p.add_argument("--power", type=float, default=1.0, help="IDW power for source-to-target remapping")
    p.add_argument("--max_dist", type=float, default=500000.0, help="Max nearest-source distance in metres")
    p.add_argument("--time_chunk", type=int, default=50)

    # Target mask
    p.add_argument("--target_mask_mode", choices=["land_only", "nominal", "none"], default="land_only")
    p.add_argument("--target_mask_var", default="status_flag")
    p.add_argument("--target_land_value", type=int, default=1)
    p.add_argument("--target_nominal_value", type=int, default=0)

    # Drift extent mask
    p.add_argument("--drift_extent_mode", choices=["nearest", "weighted_fraction", "all", "any", "none"], default="nearest")
    p.add_argument("--min_valid_weight", type=float, default=0.5)

    # Operational ice-domain drift completion
    p.add_argument(
        "--complete_drift_on_ice",
        action="store_true",
        help=(
            "Fill missing target-grid drift where non-land siconc meets the "
            "operational threshold, using a bounded nearest valid donor."
        ),
    )
    p.add_argument(
        "--ice_threshold_percent",
        type=float,
        default=15.0,
        help="Operational sea-ice concentration threshold in percent.",
    )
    p.add_argument(
        "--drift_fill_max_distance",
        type=float,
        default=300000.0,
        help="Maximum target-grid drift completion distance in metres.",
    )
    p.add_argument(
        "--drift_fill_donor_domain",
        choices=["non_land", "ice"],
        default="non_land",
        help="Use any valid non-land drift donor, or only donors over operational ice.",
    )
    p.add_argument(
        "--allow_unresolved_ice_drift",
        action="store_true",
        help=(
            "Write output even when bounded completion leaves operational ice "
            "cells without drift. The default is to fail safely."
        ),
    )
    p.add_argument(
        "--unresolved_ice_drift_policy",
        choices=["error", "exclude"],
        default="error",
        help=(
            "Fail on unresolved operational ice, or preserve model siconc "
            "while masking those cells from siconc_for_age."
        ),
    )

    # Known pole/origin repair
    p.add_argument("--fix_index_block", action="store_true", help="Discard and refill a y/x index block in siu/siv only after remapping")
    p.add_argument("--index_block_start", type=int, default=214)
    p.add_argument("--index_block_end", type=int, default=218)
    p.add_argument("--index_block_pad", type=int, default=3)
    p.add_argument("--block_method", choices=["idw", "linear"], default="idw")
    p.add_argument("--block_idw_power", type=float, default=2.0)
    p.add_argument("--no_block_nearest_fallback", action="store_true")

    # siconc
    p.add_argument("--auto_clip_siconc", action="store_true")
    p.add_argument("--clip_siconc_min", type=float, default=None)
    p.add_argument("--clip_siconc_max", type=float, default=None)

    return p.parse_args(argv)


def main(argv: list[str]) -> None:
    args = parse_args(argv)

    print("Reading source grids and metadata")
    with Dataset(args.siu, "r") as siu_ds, Dataset(args.siv, "r") as siv_ds, Dataset(args.siconc, "r") as sic_ds:
        for varname, ds, path in [
            ("siu", siu_ds, args.siu),
            ("siv", siv_ds, args.siv),
            (args.siconc_var, sic_ds, args.siconc),
        ]:
            if varname not in ds.variables:
                raise KeyError(f"{path} does not contain variable '{varname}'")

        lon_drift = nc_to_array(siu_ds.variables["longitude"])
        lat_drift = nc_to_array(siu_ds.variables["latitude"])
        lon_siv = nc_to_array(siv_ds.variables["longitude"])
        lat_siv = nc_to_array(siv_ds.variables["latitude"])
        lon_sic = nc_to_array(sic_ds.variables["longitude"])
        lat_sic = nc_to_array(sic_ds.variables["latitude"])

        dlon_siv = ((lon_drift - lon_siv + 180.0) % 360.0) - 180.0
        max_londiff_siv = np.nanmax(np.abs(dlon_siv))
        max_latdiff_siv = np.nanmax(np.abs(lat_drift - lat_siv))
        if max_londiff_siv > 1.0e-8 or max_latdiff_siv > 1.0e-8:
            print(f"Warning: siu/siv grids differ: max wrapped lon diff={max_londiff_siv}, max lat diff={max_latdiff_siv}")

        dlon_sic = ((lon_drift - lon_sic + 180.0) % 360.0) - 180.0
        max_londiff_sic = np.nanmax(np.abs(dlon_sic))
        max_latdiff_sic = np.nanmax(np.abs(lat_drift - lat_sic))
        if max_londiff_sic > 1.0e-8 or max_latdiff_sic > 1.0e-8:
            print(
                "Warning: siu/siconc lon-lat grids differ. Separate siconc weights will be used. "
                f"max wrapped lon diff={max_londiff_sic}, max lat diff={max_latdiff_sic}"
            )

        siu_var = siu_ds.variables["siu"]
        siv_var = siv_ds.variables["siv"]
        sic_var = sic_ds.variables[args.siconc_var]

        if siu_var.shape != siv_var.shape:
            raise ValueError(f"siu/siv shapes differ: {siu_var.shape} vs {siv_var.shape}")
        if siu_var.ndim != 3 or sic_var.ndim != 3:
            raise ValueError("Expected siu/siv/siconc with dimensions (time,i,j).")

        nt, ni, nj = siu_var.shape
        nt_sic, ni_sic, nj_sic = sic_var.shape

        if nt != nt_sic:
            raise ValueError(f"drift time length {nt} differs from siconc time length {nt_sic}")

    print(f"Drift source grid: {ni} x {nj}; time steps: {nt}")

    print("Computing perturbation/Jacobian native-grid basis")
    bi_e, bi_n, bj_e, bj_n, n_used, rms = compute_perturbation_basis(
        lon_drift,
        lat_drift,
        stencil_radius=args.stencil_radius,
        max_neighbor_dist=args.max_neighbor_dist,
        min_points=args.min_points,
        periodic_j=(not args.no_periodic_j),
        robust=(not args.no_robust),
    )

    valid_basis = np.isfinite(bi_e) & np.isfinite(bi_n) & np.isfinite(bj_e) & np.isfinite(bj_n)
    print(f"Valid perturbation basis cells before fill: {np.count_nonzero(valid_basis)} / {valid_basis.size}")

    if not args.no_fill_basis:
        bi_e, bi_n, bj_e, bj_n, n_filled_basis = fill_basis_nearest(bi_e, bi_n, bj_e, bj_n)
        valid_basis = np.isfinite(bi_e) & np.isfinite(bi_n) & np.isfinite(bj_e) & np.isfinite(bj_n)
        print(f"Filled basis cells: {n_filled_basis}")
        print(f"Valid perturbation basis cells after fill: {np.count_nonzero(valid_basis)} / {valid_basis.size}")

    if args.write_diagnostics:
        print(f"Writing rotation diagnostics to {args.write_diagnostics}")
        write_rotation_diagnostics(args.write_diagnostics, lon_drift, lat_drift, bi_e, bi_n, bj_e, bj_n, n_used, rms)

    print("Reading target grid")
    with Dataset(args.grid, "r") as grid_ds:
        x_tgt = nc_to_array(grid_ds.variables["x"])
        y_tgt = nc_to_array(grid_ds.variables["y"])
        lon_tgt = nc_to_array(grid_ds.variables["longitude"])
        lat_tgt = nc_to_array(grid_ds.variables["latitude"])
        target_crs = make_target_crs(grid_ds)
        target_valid = get_target_valid_mask(
            grid_ds,
            lon_tgt,
            lat_tgt,
            mode=args.target_mask_mode,
            mask_var=args.target_mask_var,
            land_value=args.target_land_value,
            nominal_value=args.target_nominal_value,
        )

    ny, nx = lon_tgt.shape
    print(f"Target grid: {ny} x {nx}")
    print(f"Target valid cells: {np.count_nonzero(target_valid)} / {target_valid.size}")

    xx_tgt, yy_tgt = np.meshgrid(x_tgt, y_tgt)

    index_inner = index_ring = None
    if args.fix_index_block:
        index_inner, index_ring, inner_bounds, window_bounds = build_index_block_masks(
            ny, nx, args.index_block_start, args.index_block_end, args.index_block_pad
        )
        index_inner &= target_valid
        index_ring &= target_valid
        print(
            f"Index-block Cartesian vector repair enabled: y/x={args.index_block_start}..{args.index_block_end} inclusive; "
            f"inner cells={np.count_nonzero(index_inner)}, ring cells={np.count_nonzero(index_ring)}, "
            f"pad={args.index_block_pad}, method={args.block_method}, idw_power={args.block_idw_power}"
        )
        print(f"  inner Python slices y={inner_bounds[0]}:{inner_bounds[1]}, x={inner_bounds[2]}:{inner_bounds[3]}")
        print(f"  window Python slices y={window_bounds[0]}:{window_bounds[1]}, x={window_bounds[2]}:{window_bounds[3]}")

    print("Projecting source grids to target CRS")
    transformer = Transformer.from_crs(CRS.from_epsg(4326), target_crs, always_xy=True)
    x_src_drift, y_src_drift = transformer.transform(lon_drift, lat_drift)
    x_src_sic, y_src_sic = transformer.transform(lon_sic, lat_sic)

    x_src_drift_flat = np.asarray(x_src_drift, dtype=np.float64).ravel()
    y_src_drift_flat = np.asarray(y_src_drift, dtype=np.float64).ravel()
    x_src_sic_flat = np.asarray(x_src_sic, dtype=np.float64).ravel()
    y_src_sic_flat = np.asarray(y_src_sic, dtype=np.float64).ravel()

    x_tgt_flat = xx_tgt.ravel()
    y_tgt_flat = yy_tgt.ravel()
    target_valid_flat = target_valid.ravel()

    drift_source_good = np.isfinite(x_src_drift_flat) & np.isfinite(y_src_drift_flat)
    sic_source_good = np.isfinite(x_src_sic_flat) & np.isfinite(y_src_sic_flat)

    if args.weights and os.path.exists(args.weights):
        print(f"Loading weights from {args.weights}")
        with np.load(args.weights) as npz:
            required = {"drift_indices", "drift_weights", "drift_nearest_index", "sic_indices", "sic_weights", "sic_nearest_index"}
            if not required.issubset(set(npz.files)):
                raise ValueError("Weights file is not compatible. Use a new --weights filename or delete the old file.")

            drift_indices = npz["drift_indices"]
            drift_weights = npz["drift_weights"]
            drift_nearest_index = npz["drift_nearest_index"]
            sic_indices = npz["sic_indices"]
            sic_weights = npz["sic_weights"]
            sic_nearest_index = npz["sic_nearest_index"]

        if drift_indices.shape[0] != target_valid_flat.size or sic_indices.shape[0] != target_valid_flat.size:
            raise ValueError("Weight file target size does not match current target grid.")
    else:
        print("Computing drift source-to-target weights")
        drift_indices, drift_weights, drift_nearest_index = build_weights(
            x_src_drift_flat,
            y_src_drift_flat,
            drift_source_good,
            x_tgt_flat,
            y_tgt_flat,
            target_valid_flat,
            k=args.k,
            max_dist=args.max_dist,
            power=args.power,
        )

        print("Computing siconc source-to-target weights")
        sic_indices, sic_weights, sic_nearest_index = build_weights(
            x_src_sic_flat,
            y_src_sic_flat,
            sic_source_good,
            x_tgt_flat,
            y_tgt_flat,
            target_valid_flat,
            k=args.k,
            max_dist=args.max_dist,
            power=args.power,
        )

        if args.weights:
            np.savez(
                args.weights,
                drift_indices=drift_indices,
                drift_weights=drift_weights,
                drift_nearest_index=drift_nearest_index,
                sic_indices=sic_indices,
                sic_weights=sic_weights,
                sic_nearest_index=sic_nearest_index,
                drift_source_shape=np.asarray(lon_drift.shape, dtype=np.int64),
                sic_source_shape=np.asarray(lon_sic.shape, dtype=np.int64),
                target_shape=np.asarray(lon_tgt.shape, dtype=np.int64),
                k=np.asarray(args.k, dtype=np.int64),
                power=np.asarray(args.power, dtype=np.float64),
                max_dist=np.asarray(args.max_dist, dtype=np.float64),
                target_mask_mode=np.asarray(args.target_mask_mode),
            )
            print(f"Saved weights to {args.weights}")

    if args.auto_clip_siconc:
        sic_clip_min = 0.0 if args.clip_siconc_min is None else args.clip_siconc_min
        sic_clip_max = 100.0 if args.clip_siconc_max is None else args.clip_siconc_max
    else:
        sic_clip_min = args.clip_siconc_min
        sic_clip_max = args.clip_siconc_max

    print(f"Creating output {args.out}")
    out = create_output(
        args.out,
        args.siu,
        args.siv,
        args.siconc,
        args.grid,
        nt,
        ny,
        nx,
        args.siconc_var,
        args.target_mask_mode,
    )

    out_siu = out.variables["siu"]
    out_siv = out.variables["siv"]
    out_sic = out.variables[args.siconc_var]
    out_age_sic = out.variables["siconc_for_age"]
    out_drift_support = out.variables["drift_supported"]

    total_block_fixed = 0
    total_ice_cells = 0
    total_missing_ice_drift = 0
    total_drift_completed = 0
    total_unresolved_ice_drift = 0
    total_excluded_ice_drift = 0
    drift_completion_distances_m = []
    minimum_unresolved_distance_m = np.nan

    print("Processing time steps")
    with Dataset(args.siu, "r") as siu_ds, Dataset(args.siv, "r") as siv_ds, Dataset(args.siconc, "r") as sic_ds:
        siu_in = siu_ds.variables["siu"]
        siv_in = siv_ds.variables["siv"]
        sic_in = sic_ds.variables[args.siconc_var]

        for t0 in range(0, nt, args.time_chunk):
            t1 = min(nt, t0 + args.time_chunk)
            count = t1 - t0
            print(f"  time {t0}..{t1 - 1} / {nt - 1}")

            siu_block_native = read_time_block_as_float(siu_in, t0, count)
            siv_block_native = read_time_block_as_float(siv_in, t0, count)
            sic_block_native = read_time_block_as_float(sic_in, t0, count)

            siu_block_out = np.full((count, ny, nx), np.nan, dtype=np.float32)
            siv_block_out = np.full((count, ny, nx), np.nan, dtype=np.float32)
            sic_block_out = np.full((count, ny, nx), np.nan, dtype=np.float32)
            age_sic_block_out = np.full(
                (count, ny, nx), np.nan, dtype=np.float32
            )
            drift_support_block_out = np.full(
                (count, ny, nx), np.uint8(255), dtype=np.uint8
            )

            for kk in range(count):
                u_east_native, v_north_native = rotate_native_to_east_north(
                    siu_block_native[kk],
                    siv_block_native[kk],
                    bi_e,
                    bi_n,
                    bj_e,
                    bj_n,
                )

                # Drift: true vector remapping via Cartesian tangent vectors.
                siu_laea, siv_laea = vector_cartesian_remap(
                    u_east_native,
                    v_north_native,
                    lon_drift,
                    lat_drift,
                    lon_tgt,
                    lat_tgt,
                    drift_indices,
                    drift_weights,
                )

                # Scalar concentration remapping.
                sic_laea = apply_weights(sic_block_native[kk].ravel(), sic_indices, sic_weights).reshape(ny, nx)

                # Drift extent mask prevents extrapolation beyond source finite-drift area.
                drift_extent_flat = make_drift_extent_mask(
                    u_east_native.ravel(),
                    v_north_native.ravel(),
                    drift_indices,
                    drift_weights,
                    drift_nearest_index,
                    mode=args.drift_extent_mode,
                    min_valid_weight=args.min_valid_weight,
                )
                drift_extent = drift_extent_flat.reshape(ny, nx) & target_valid

                siu_laea[~drift_extent] = np.nan
                siv_laea[~drift_extent] = np.nan
                sic_laea[~target_valid] = np.nan

                if sic_clip_min is not None or sic_clip_max is not None:
                    sic_laea = np.clip(
                        sic_laea,
                        -np.inf if sic_clip_min is None else sic_clip_min,
                        np.inf if sic_clip_max is None else sic_clip_max,
                    )

                # Optional final repair for the known pole/origin index block.
                if args.fix_index_block:
                    siu_laea, siv_laea, nfix, _, _ = cartesian_vector_block_repair(
                        siu_laea,
                        siv_laea,
                        xx_tgt,
                        yy_tgt,
                        lon_tgt,
                        lat_tgt,
                        index_inner,
                        index_ring,
                        method=args.block_method,
                        idw_power=args.block_idw_power,
                        nearest_fallback=(not args.no_block_nearest_fallback),
                    )
                    total_block_fixed += nfix

                if args.complete_drift_on_ice:
                    siu_laea, siv_laea, completion = complete_drift_where_ice(
                        siu_laea,
                        siv_laea,
                        sic_laea,
                        xx_tgt,
                        yy_tgt,
                        lon_tgt,
                        lat_tgt,
                        target_valid,
                        concentration_threshold_percent=args.ice_threshold_percent,
                        maximum_distance_m=args.drift_fill_max_distance,
                        donor_domain=args.drift_fill_donor_domain,
                    )
                    total_ice_cells += completion["ice_cells"]
                    total_missing_ice_drift += completion["missing_cells"]
                    total_drift_completed += completion["repaired_cells"]
                    total_unresolved_ice_drift += completion["unresolved_cells"]
                    if completion["repair_distances_m"].size:
                        drift_completion_distances_m.append(
                            completion["repair_distances_m"]
                        )
                    unresolved_distance = completion[
                        "minimum_unresolved_distance_m"
                    ]
                    if np.isfinite(unresolved_distance):
                        minimum_unresolved_distance_m = (
                            unresolved_distance
                            if not np.isfinite(minimum_unresolved_distance_m)
                            else min(
                                minimum_unresolved_distance_m,
                                unresolved_distance,
                            )
                        )
                    if completion["missing_cells"]:
                        absolute_time_index = t0 + kk
                        distances = completion["repair_distances_m"]
                        maximum_repair_km = (
                            float(np.max(distances)) / 1000.0
                            if distances.size else 0.0
                        )
                        minimum_unresolved_km = (
                            float(unresolved_distance) / 1000.0
                            if np.isfinite(unresolved_distance) else np.nan
                        )
                        print(
                            "    drift completion "
                            f"time_index={absolute_time_index}: "
                            f"ice_cells={completion['ice_cells']}, "
                            f"missing={completion['missing_cells']}, "
                            f"repaired={completion['repaired_cells']}, "
                            f"unresolved={completion['unresolved_cells']}, "
                            f"maximum_repair_distance={maximum_repair_km:.2f} km, "
                            "minimum_unresolved_distance="
                            f"{minimum_unresolved_km:.2f} km"
                        )
                    if (
                        completion["unresolved_cells"]
                        and args.unresolved_ice_drift_policy == "error"
                        and not args.allow_unresolved_ice_drift
                    ):
                        out.drift_completion_enabled = np.int8(1)
                        out.drift_completion_complete = np.int8(0)
                        out.drift_completion_threshold_percent = np.float64(
                            args.ice_threshold_percent
                        )
                        out.drift_completion_max_distance_km = np.float64(
                            args.drift_fill_max_distance / 1000.0
                        )
                        out.drift_completion_repaired_cell_times = np.int64(
                            total_drift_completed
                        )
                        out.drift_completion_unresolved_cell_times = np.int64(
                            total_unresolved_ice_drift
                        )
                        out.drift_completion_min_unresolved_distance_km = (
                            np.float64(minimum_unresolved_km)
                        )
                        out.sync()
                        out.close()
                        raise RuntimeError(
                            f"Time index {t0 + kk}: "
                            f"{completion['unresolved_cells']} non-land target "
                            "cells at or above the operational sea-ice "
                            "concentration threshold still lack finite siu/siv "
                            "after bounded drift completion. Partial output was "
                            "not marked complete. Nearest unresolved donor "
                            f"distance={minimum_unresolved_km:.2f} km; configured "
                            "maximum="
                            f"{args.drift_fill_max_distance / 1000.0:.2f} km."
                        )

                (
                    age_sic_laea,
                    drift_supported,
                    unresolved_ice,
                ) = make_age_domain_fields(
                    sic_laea,
                    siu_laea,
                    siv_laea,
                    target_valid,
                    args.ice_threshold_percent,
                    args.unresolved_ice_drift_policy,
                )
                if args.unresolved_ice_drift_policy == "exclude":
                    total_excluded_ice_drift += int(
                        np.count_nonzero(unresolved_ice)
                    )
                siu_block_out[kk, :, :] = np.where(np.isfinite(siu_laea), siu_laea, FILL_VALUE).astype(np.float32)
                siv_block_out[kk, :, :] = np.where(np.isfinite(siv_laea), siv_laea, FILL_VALUE).astype(np.float32)
                sic_block_out[kk, :, :] = np.where(np.isfinite(sic_laea), sic_laea, FILL_VALUE).astype(np.float32)
                age_sic_block_out[kk, :, :] = np.where(
                    np.isfinite(age_sic_laea),
                    age_sic_laea,
                    FILL_VALUE,
                ).astype(np.float32)
                drift_support_block_out[kk, target_valid] = (
                    drift_supported[target_valid].astype(np.uint8)
                )

            out_siu[t0:t1, :, :] = siu_block_out
            out_siv[t0:t1, :, :] = siv_block_out
            out_sic[t0:t1, :, :] = sic_block_out
            out_age_sic[t0:t1, :, :] = age_sic_block_out
            out_drift_support[t0:t1, :, :] = drift_support_block_out

    out.remapping_settings = (
        f"source-to-target IDW weights: k={args.k}, power={args.power}, max_dist={args.max_dist} m; "
        "drift remapped as Cartesian tangent vectors; siconc remapped as scalar."
    )
    out.perturbation_rotation_settings = (
        f"stencil_radius={args.stencil_radius}, max_neighbor_dist={args.max_neighbor_dist} m, "
        f"min_points={args.min_points}, robust={not args.no_robust}, periodic_j={not args.no_periodic_j}."
    )
    out.drift_extent_method = (
        f"mode={args.drift_extent_mode}, min_valid_weight={args.min_valid_weight}; "
        "applied only to siu/siv."
    )
    out.index_block_repair = (
        f"enabled={args.fix_index_block}; y/x={args.index_block_start}..{args.index_block_end}; "
        f"pad={args.index_block_pad}; method={args.block_method}; idw_power={args.block_idw_power}; "
        f"total_fixed_cell_times={total_block_fixed}; applied to siu/siv only after remapping."
    )
    completion_distances = (
        np.concatenate(drift_completion_distances_m)
        if drift_completion_distances_m
        else np.empty(0, dtype=np.float64)
    )
    out.drift_completion_enabled = np.int8(args.complete_drift_on_ice)
    out.drift_completion_complete = np.int8(
        args.complete_drift_on_ice
        and (
            total_unresolved_ice_drift == 0
            or args.unresolved_ice_drift_policy == "exclude"
        )
    )
    out.drift_completion_unresolved_policy = args.unresolved_ice_drift_policy
    out.drift_completion_threshold_percent = np.float64(
        args.ice_threshold_percent
    )
    out.drift_completion_max_distance_km = np.float64(
        args.drift_fill_max_distance / 1000.0
    )
    out.drift_completion_donor_domain = args.drift_fill_donor_domain
    out.drift_completion_ice_cell_times = np.int64(total_ice_cells)
    out.drift_completion_missing_cell_times = np.int64(
        total_missing_ice_drift
    )
    out.drift_completion_repaired_cell_times = np.int64(
        total_drift_completed
    )
    out.drift_completion_unresolved_cell_times = np.int64(
        total_unresolved_ice_drift
    )
    out.drift_completion_excluded_cell_times = np.int64(
        total_excluded_ice_drift
    )
    out.drift_exclusion_complete = np.int8(
        args.unresolved_ice_drift_policy == "exclude"
        and total_excluded_ice_drift == total_unresolved_ice_drift
    )
    out.drift_completion_max_repair_distance_km = np.float64(
        np.max(completion_distances) / 1000.0
        if completion_distances.size else 0.0
    )
    out.drift_completion_mean_repair_distance_km = np.float64(
        np.mean(completion_distances) / 1000.0
        if completion_distances.size else 0.0
    )
    out.drift_completion_p95_repair_distance_km = np.float64(
        np.percentile(completion_distances, 95.0) / 1000.0
        if completion_distances.size else 0.0
    )
    out.drift_completion_min_unresolved_distance_km = np.float64(
        minimum_unresolved_distance_m / 1000.0
        if np.isfinite(minimum_unresolved_distance_m) else np.nan
    )
    out.drift_completion_method = (
        "Missing drift on non-land target cells at or above the operational "
        "siconc threshold was completed from the nearest bounded target-grid "
        f"donor in domain={args.drift_fill_donor_domain}. Donor east/north "
        "vectors were converted to Earth-centred Cartesian tangent vectors "
        "and projected onto each recipient's local basis. Unresolved cells "
        f"used policy={args.unresolved_ice_drift_policy}; exclusion preserves "
        "model siconc and masks only siconc_for_age."
    )
    print(
        "Drift completion summary: "
        f"enabled={args.complete_drift_on_ice}, "
        f"ice_cell_times={total_ice_cells}, "
        f"missing={total_missing_ice_drift}, "
        f"repaired={total_drift_completed}, "
        f"unresolved={total_unresolved_ice_drift}, "
        f"excluded_from_age={total_excluded_ice_drift}, "
        "maximum_repair_distance="
        f"{out.drift_completion_max_repair_distance_km:.2f} km, "
        "p95_repair_distance="
        f"{out.drift_completion_p95_repair_distance_km:.2f} km"
    )

    out.close()
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main(sys.argv[1:])
