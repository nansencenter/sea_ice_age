import glob
from datetime import timedelta
from multiprocessing import Pool
import os

from matplotlib.tri import Triangulation
import numpy as np
from shapely.geometry import Polygon
from scipy.interpolate import RectBivariateSpline
from scipy.spatial import KDTree
from scipy.ndimage import distance_transform_edt

from .remeshing import (
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

# TODO:
# Outdated function - remove later
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

class GetWeights:
    def __init__(self, xa, ya, ta, xo, yo, to, treea, search_dist):
        self.xa = xa
        self.ya = ya
        self.ta = ta
        self.xo = xo
        self.yo = yo
        self.to = to
        self.treea = treea
        self.search_dist = search_dist

    def __call__(self, i):
        src2dst = []
        weights = []
        x = self.xo[self.to[i]].mean()
        y = self.yo[self.to[i]].mean()
        near_elems_1 = self.treea.query_ball_point([x, y], self.search_dist)
        po = make_polygon_shapely(self.xo[self.to[i]], self.yo[self.to[i]])
        for j in near_elems_1:
            pa = make_polygon_shapely(self.xa[self.ta[j]], self.ya[self.ta[j]])
            p_int = po.intersection(pa)
            if p_int.area > 0:
                src2dst.append((j, i))
                weights.append(p_int.area / po.area)
        return src2dst, weights

def compute_mapping(tri_a, tri_o, search_dist, cores=4):
    xa, ya, ta = tri_a.x, tri_a.y, tri_a.triangles
    xo, yo, to = tri_o.x, tri_o.y, tri_o.triangles


    same_in_a, same_in_o = get_same_elements(xa[ta], xo[to])
    new_elems = np.ones(tri_o.triangles.shape[0], dtype=bool)
    new_elems[same_in_o] = False
    new_elem_idx = np.nonzero(new_elems)[0]

    src2dst = np.column_stack([same_in_a, same_in_o]).tolist()
    weights = [1] * len(src2dst)
    elem_xa = xa[ta].mean(axis=1)
    elem_ya = ya[ta].mean(axis=1)
    get_src2dst_weights = GetWeights(xa, ya, ta, xo, yo, to,
                                     KDTree(np.column_stack([elem_xa, elem_ya])),
                                     search_dist)

    with Pool(cores) as p:
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

class InterpolateSic:
    def __init__(self, sid_dir, xc, yc, landmask, min_lm, force):
        self.sid_dir = sid_dir
        self.xc = xc
        self.yc = yc
        self.landmask = landmask
        self.min_lm = min_lm
        self.force = force

    #for mesh_file in mesh_files:
    def __call__(self, mesh_file):
        basename = os.path.basename(mesh_file)
        sic_dst_file = mesh_file.replace("mesh_", "sic_").replace("/mesh/", "/sic/")
        if os.path.exists(sic_dst_file) and not self.force:
            return
        print(sic_dst_file)
        mesh_date = basename.split('.')[0].split('_')[-1]
        mesh_year = mesh_date[:4]
        # load mesh
        with np.load(mesh_file) as d:
            x = d['x']
            y = d['y']
            t = d['t']

        # load SIC
        file_mask = f'{self.sid_dir}/{mesh_year}/ice_drift_nh_ease*{mesh_date}1200.nc.npz'
        try:
            sic_src_file = glob.glob(file_mask)[0]
        except:
            raise ValueError(f'Cannot find {file_mask}')
        try:
            cgrd = np.load(sic_src_file)['c']
        except:
            print(f'Cannot load c from {sic_src_file}')
            raise ValueError
        cgrd[np.isnan(cgrd)] = 0

        # interpolate SIC
        try:
            c = RectBivariateSpline(self.xc, self.yc, cgrd[::-1], kx=1, ky=1)(y[t].mean(axis=1), x[t].mean(axis=1), grid=False)
        except:
            raise ValueError(f'Fail to interpolate SIC in {mesh_file}')
        # interpolate landmask
        try:
            lm = RectBivariateSpline(self.xc, self.yc, self.landmask[::-1], kx=1, ky=1)(y[t].mean(axis=1), x[t].mean(axis=1), grid=False)
        except:
            raise ValueError(f'Fail to interpolate LM in {mesh_file}')
        # clear land elements
        c[lm > self.min_lm] = 0
        os.makedirs(os.path.split(sic_dst_file)[0], exist_ok=True)
        np.savez(sic_dst_file, c=c)


class IrregularGridInterpolator(object):
    def __init__(self, x0, y0, x1, y1, triangles=None):
        '''
        Parameters:
        -----------
        x0 : np.ndarray(float)
            x-coords of source points
        y0 : np.ndarray(float)
            y-coords of source points
        x1 : np.ndarray(float)
            x-coords of destination points
        y1 : np.ndarray(float)
            y-coords of destination points
        triangles : np.ndarray(int)
            shape (num_triangles, 3)
            indices of nodes for each triangle

        Sets:
        -----
        self.inside: np.ndarray(bool)
            shape = (num_target_points,)
        self.vertices: np.ndarray(int)
            shape = (num_good_target_points, 3)
            good target points are those inside the source triangulation
        self.weights: np.ndarray(float)
            shape = (num_good_target_points, 3)
            good target points are those inside the source triangulation

        Follows this suggestion:
        https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids
        x_target[i] = sum_{j=0}^2 weights[i, j]*x_source[vertices[i, j]]
        y_target[i] = sum_{j=0}^2 weights[i, j]*y_source[vertices[i, j]]
        We can do (linear) interpolation by replacing x_target, x_source with z_target, z_source
        where z_source is the field to be interpolated and z_target is the interpolated field
        '''

        # define and triangulate source points
        self.src_shape = x0.shape
        self.src_points = np.array([x0.flatten(), y0.flatten()]).T
        self.tri = Triangulation(x0.flatten(), y0.flatten(), triangles=triangles)
        self.tri_finder = self.tri.get_trifinder()
        self.num_triangles = len(self.tri.triangles)
        self._set_transform()

        # define target points
        self.dst_points = np.array([x1.flatten(), y1.flatten()]).T
        self.dst_shape = x1.shape
        self.triangle_map = self.tri_finder(x1, y1)
        self.dst_mask = (self.triangle_map < 0)
        self.triangle_map[self.dst_mask] = 0
        self.inside = ~self.dst_mask.flatten()

        """
        get barycentric coords
        https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Barycentric_coordinates_on_triangles
        each row of bary is (lambda_1, lambda_2) for 1 destination point
        """
        d = 2
        inds = self.triangle_map.flatten()[self.inside]
        self.vertices = np.take(self.tri.triangles, inds, axis=0)
        temp = np.take(self.transform, inds, axis=0)
        delta = self.dst_points[self.inside] - temp[:, d]
        bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)

        # set weights
        self.weights = np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

    def _set_transform(self):
        """
        Used for getting the barycentric coordinates on a triangle.
        Follows:
        https://en.wikipedia.org/wiki/Barycentric_coordinate_system#Barycentric_coordinates_on_triangles

        Sets:
        -----
        self.transform : numpy.ndarray
            For the i-th triangle,
                self.transform[i] = [[a', b'], [c', d'], [x_3, y3]]
            where the first 2 rows are the inverse of the matrix T in the wikipedia link
            and (x_3, y_3) are the coordinates of the 3rd vertex of the triangle
        """
        x = self.tri.x[self.tri.triangles]
        y = self.tri.y[self.tri.triangles]
        a = x[:,0] - x[:,2]
        b = x[:,1] - x[:,2]
        c = y[:,0] - y[:,2]
        d = y[:,1] - y[:,2]
        det = a*d-b*c

        self.transform = np.zeros((self.num_triangles, 3, 2))
        self.transform[:,0,0] = d/det
        self.transform[:,0,1] = -b/det
        self.transform[:,1,0] = -c/det
        self.transform[:,1,1] = a/det
        self.transform[:,2,0] = x[:,2]
        self.transform[:,2,1] = y[:,2]

    def interp_field(self, fld, method='linear'):
        """
        Interpolate field from elements elements or nodes of source triangulation
        to destination points

        Parameters:
        -----------
        fld: np.ndarray
            field to be interpolated
        method : str
            interpolation method if interpolating from nodes
            - 'linear'  : linear interpolation
            - 'nearest' : nearest neighbour

        Returns:
        -----------
        fld_interp : np.ndarray
            field interpolated onto the destination points
        """
        if fld.shape == self.src_shape:
            return self._interp_nodes(fld, method=method)
        fld_ = fld.flatten()
        if len(fld_) == self.num_triangles:
            return self._interp_elements(fld_)
        msg = f"""Field to interpolate should have the same size as the source points
        i.e. {self.src_shape}, or be a vector with the same number of triangles
        as the source triangulation i.e. self.num_triangles"""
        raise ValueError(msg)

    def _interp_elements(self, fld):
        """
        Interpolate field from elements of source triangulation to destination points

        Parameters:
        -----------
        fld: np.ndarray
            field to be interpolated

        Returns:
        -----------
        fld_interp : np.ndarray
            field interpolated onto the destination points
        """
        fld_interp = fld[self.triangle_map]
        fld_interp[self.dst_mask] = np.nan
        return fld_interp

    def _interp_nodes(self, fld, method='linear'):
        """
        Interpolate field from nodes of source triangulation to destination points

        Parameters:
        -----------
        fld: np.ndarray
            field to be interpolated
        method : str
            interpolation method
            - 'linear'  : linear interpolation
            - 'nearest' : nearest neighbour

        Returns:
        -----------
        fld_interp : np.ndarray
            field interpolated onto the destination points
        """
        ndst = self.dst_points.shape[0]
        fld_interp = np.full((ndst,), np.nan)
        w = self.weights
        if method == 'linear':
            # sum over the weights for each node of triangle
            v = self.vertices # shape = (ngood,3)
            fld_interp[self.inside] = np.einsum(
                    'nj,nj->n', np.take(fld.flatten(), v), w)

        elif method == 'nearest':
            # find the node of the triangle with the maximum weight
            v = np.array(self.vertices) # shape = (ngood,3)
            v = v[np.arange(len(w), dtype=int), np.argmax(w, axis=1)] # shape = (ngood,)
            fld_interp[self.inside] = fld.flatten()[v]

        else:
            raise ValueError("'method' should be 'nearest' or 'linear'")

        return fld_interp.reshape(self.dst_shape)
