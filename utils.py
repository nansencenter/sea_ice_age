from multiprocessing import Pool
import os
from datetime import timedelta

from matplotlib.tri import Triangulation
import numpy as np
from shapely.geometry import Polygon
from scipy.interpolate import RectBivariateSpline
from scipy.spatial import KDTree
from scipy.ndimage import distance_transform_edt

from remeshing import (
    get_area,
    remove_small_elements,
    collapse_short_edges,
    remove_hanging_elements,
    collapse_small_angle_edges,
    remove_long_edges,
    get_same_elements,
    measure,
    laplace_smooth,
    clean_mesh,
    )

def make_polygon_shapely(x, y):
    return Polygon([(x[0],y[0]), (x[1],y[1]), (x[2],y[2])])

def resample_velocities(u, v, m, xc, yc, x, y, max_mask_fix_node, max_mask_zero_speed):
    u0 = RectBivariateSpline(xc[1:-1:3], yc[1:-1:3], u[::-1], kx=1, ky=1)(y, x, grid=False)
    v0 = RectBivariateSpline(xc[1:-1:3], yc[1:-1:3], v[::-1], kx=1, ky=1)(y, x, grid=False)
    mask = RectBivariateSpline(xc, yc, m[::-1], kx=1, ky=1)(y, x, grid=False)
    fixed_nodes_idx = np.nonzero(mask < max_mask_fix_node)[0]
    zero_speed = np.nonzero(mask < max_mask_zero_speed)[0]
    u0[zero_speed] = 0
    v0[zero_speed] = 0
    return u0, v0, fixed_nodes_idx

def get_distance_weights(x0, y0, t0, a1, max_dist=50):
    neg_elem_x = x0[t0[a1<0]].mean(axis=1)
    neg_elem_y = y0[t0[a1<0]].mean(axis=1)
    neg_trr = KDTree(np.vstack([neg_elem_x, neg_elem_y]).T)
    dist, _ = neg_trr.query(np.vstack([x0, y0]).T, k=1)
    weights = 1 - dist / max_dist
    weights[weights < 0] = 0
    return weights

def decrease_speed_factor(weights, max_factor):
    return (1 - weights * max_factor)

def advect_nodes(tri0, u0, v0, max_dist0):
    x0, y0, t0 = tri0.x, tri0.y, tri0.triangles
    x1 = x0 + u0
    y1 = y0 + v0
    a1 = get_area(x1, y1, t0)
    if a1.min() < 0:
        for max_dist in [max_dist0, max_dist0*2, max_dist0*3, max_dist0*4]:
            for max_factor in [0.25, 0.5, 0.75, 0.9]:
                weights = get_distance_weights(x0, y0, t0, a1, max_dist)
                factor = decrease_speed_factor(weights, max_factor)
                u0a = u0 * factor
                v0a = v0 * factor
                x1a = x0 + u0a
                y1a = y0 + v0a
                a1a = get_area(x1a, y1a, t0)
                x1, y1, a1 = x1a, y1a, a1a
                if a1a.min() > 0:
                    print('Fixed by negative area factor', max_factor, max_dist, factor[factor != 1].size)
                    break
            if a1a.min() > 0:
                break
    if a1.min() < 0:
        print('NEGATIVE AREA!')
    return Triangulation(x1, y1, tri0.triangles)

def get_mesh_files(idate, lag_dir, mesh_init_file):
    mesh_src_dir = idate.strftime(f'{lag_dir}/%Y')
    mesh_src_file = idate.strftime(f'{mesh_src_dir}/mesh_%Y%m%d.npz')
    if not os.path.exists(mesh_src_file):
        print('Load INIT')
        mesh_src_file = mesh_init_file

    mesh_dst_date = idate + timedelta(1)
    mesh_dst_dir = mesh_dst_date.strftime(f'{lag_dir}/%Y')
    os.makedirs(mesh_dst_dir, exist_ok=True)
    mesh_dst_file = mesh_dst_date.strftime(f'{mesh_dst_dir}/mesh_%Y%m%d.npz')

    return mesh_src_file, mesh_dst_file

def load_data(ifile, mesh_src_file):
    with np.load(ifile) as d:
        u = d['u']
        v = -d['v']
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0
    with np.load(mesh_src_file) as d:
        x0 = d['x']
        y0 = d['y']
        t0 = d['t']
    return u, v, Triangulation(x0, y0, t0)

def remesh(tri_a, fixed_nodes_idx, min_area, min_edge_length, min_edge_angle, max_edge_length, verbose=False):
    tri, n = remove_small_elements(tri_a, fixed_nodes_idx, min_area)
    if verbose: print(f'Removed {n} small elements')

    tri, n = collapse_short_edges(tri, min_edge_length, fixed_nodes_idx=fixed_nodes_idx)
    tri = remove_hanging_elements(tri)
    if verbose: print(f'Collapsed {n} short edges')

    tri, n = collapse_small_angle_edges(tri, min_edge_angle, fixed_nodes_idx=fixed_nodes_idx)
    tri = remove_hanging_elements(tri)
    if verbose: print(f'Collapsed {n} small angle edges')

    tri, n = remove_long_edges(tri, max_edge_length, fixed_nodes_idx=fixed_nodes_idx)
    if verbose: print(f'Removed {n} long edges')
    return tri

def optimize_mesh(tri_a, tri_r, fixed_nodes_idx, min_ap_ratio=0.18, min_area=20, verbose=False):
    _, same_in_r = get_same_elements(tri_a.triangles, tri_r.triangles)
    new_elems = np.ones(tri_r.triangles.shape[0], dtype=bool)
    new_elems[same_in_r] = False
    _, _, _, ap_r = measure(tri_r.x, tri_r.y, tri_r.triangles)
    new_elems[ap_r < min_ap_ratio] = True
    new_elems[tri_r.neighbors[new_elems].flatten()] = True
    new_elems[tri_r.neighbors[new_elems].flatten()] = True

    fixed_nodes = set(list(np.unique(tri_r.triangles))) - set(list(np.unique(tri_r.triangles[new_elems])))
    fixed_nodes = fixed_nodes.union(set(list(fixed_nodes_idx)))
    fixed_nodes = np.array(list(fixed_nodes))

    tri = laplace_smooth(tri_r, iterations=40, factor=0.5, fixed_nodes=fixed_nodes)

    tri, n = remove_small_elements(tri, fixed_nodes_idx, min_area=min_area)
    if verbose: print(f'Removed {n} small elements')
    tri = clean_mesh(tri)
    return tri

def get_src2dst_weights(i):
    xo = get_src2dst_weights.xo
    yo = get_src2dst_weights.yo
    to = get_src2dst_weights.to
    xa = get_src2dst_weights.xa
    ya = get_src2dst_weights.ya
    ta = get_src2dst_weights.ta
    treea = get_src2dst_weights.treea
    search_dist = get_src2dst_weights.search_dist

    src2dst = []
    weights = []
    x = xo[to[i]].mean()
    y = yo[to[i]].mean()
    near_elems_1 = treea.query_ball_point([x, y], search_dist)
    po = make_polygon_shapely(xo[to[i]], yo[to[i]])
    for j in near_elems_1:
        pa = make_polygon_shapely(xa[ta[j]], ya[ta[j]])
        p_int = po.intersection(pa)
        if p_int.area > 0:
            src2dst.append((j, i))
            weights.append(p_int.area / po.area)
    return src2dst, weights

def compute_mapping(tri_a, tri_o, search_dist):
    global get_src2dst_weights
    xa, ya, ta = tri_a.x, tri_a.y, tri_a.triangles
    xo, yo, to = tri_o.x, tri_o.y, tri_o.triangles
    get_src2dst_weights.xa = xa
    get_src2dst_weights.ya = ya
    get_src2dst_weights.ta = ta
    get_src2dst_weights.xo = xo
    get_src2dst_weights.yo = yo
    get_src2dst_weights.to = to
    get_src2dst_weights.search_dist = search_dist

    same_in_a, same_in_o = get_same_elements(xa[ta], xo[to])
    new_elems = np.ones(tri_o.triangles.shape[0], dtype=bool)
    new_elems[same_in_o] = False
    new_elem_idx = np.nonzero(new_elems)[0]

    src2dst = np.column_stack([same_in_a, same_in_o]).tolist()
    weights = [1] * len(src2dst)
    elem_xa = xa[ta].mean(axis=1)
    elem_ya = ya[ta].mean(axis=1)
    get_src2dst_weights.treea = KDTree(np.column_stack([elem_xa, elem_ya]))

    with Pool(5) as p:
        res = p.map(get_src2dst_weights, new_elem_idx.tolist())
        for r in res:
            src2dst.extend(r[0])
            weights.extend(r[1])
    src2dst = np.array(src2dst)
    weights = np.array(weights)
    return src2dst, weights

def get_area_ratio(tri0, tri_a, tri_o, src2dst, weights):
    area0 = get_area(tri0.x, tri0.y, tri0.triangles)
    area1 = get_area(tri_a.x, tri_a.y, tri_a.triangles)
    area_ratio1 = area1 / area0
    area_ratio7 = np.zeros(tri_o.triangles.shape[0]) + 0
    np.add.at(area_ratio7, src2dst[:,1], area_ratio1[src2dst[:,0]] * weights)
    return area_ratio7

def fill_gaps(array, mask, distance=15):
    """ Fill gaps in input raster

    Parameters
    ----------
    array : 2D numpy.array
        Ratser with deformation field
    mask : 2D numpy.array
        Where are gaps
    distance : int
        Minimum size of gap to fill

    Returns
    -------
    arra : 2D numpy.array
        Ratser with gaps filled

    """
    dist, indi = distance_transform_edt(
        mask,
        return_distances=True,
        return_indices=True)
    gpi = dist <= distance
    r,c = indi[:,gpi]
    new_array = np.array(array)
    new_array[gpi] = array[r,c]
    return new_array