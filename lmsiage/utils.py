import glob
from datetime import timedelta
from multiprocessing import Pool
import os

from matplotlib.tri import Triangulation
import numpy as np
from numba import njit, prange
from scipy.interpolate import RectBivariateSpline
from scipy.spatial import KDTree
from scipy.ndimage import distance_transform_edt

from .mesh_file import MeshFile

def jacobian(x0, y0, x1, y1, x2, y2):
    return (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)

def get_area(x, y, t):
    return .5*jacobian(x[t][:,0], y[t][:,0], x[t][:,1], y[t][:,1], x[t][:,2], y[t][:,2])

def measure(x, y, t):
    dx = np.diff(np.hstack([x[t], x[t][:,0][None].T]))
    dy = np.diff(np.hstack([y[t], y[t][:,0][None].T]))
    edges = np.hypot(dx, dy)
    perim = edges.sum(axis=1)
    area = get_area(x, y, t)
    ap_ratio = area**0.5/ perim
    return area, edges, perim, ap_ratio

def clean_mesh(tri):
    """
    Remove unused points from a mesh and reindex triangulation

    Parameters:
    x, y (numpy arrays) - coordinates of all points
    t (numpy array) - triangulation indices

    Returns:
    x_new, y_new, t_new - cleaned mesh data
    """
    x, y, t = tri.x, tri.y, tri.triangles
    # Find all unique indices that are used in the triangulation
    used_indices = np.unique(t)

    # Create a mapping from old indices to new indices
    index_map = np.zeros(len(x), dtype=np.int32) - 1
    index_map[used_indices] = np.arange(len(used_indices))

    # Create new coordinate arrays containing only used points
    x_new = x[used_indices]
    y_new = y[used_indices]

    # Reindex triangulation
    t_new = index_map[t]

    #print(f"Removed {len(x) - len(x_new)} unused points from mesh")
    return Triangulation(x_new, y_new, t_new)

def get_same_elements(t0, t1):
    set0 = set([tuple(i) for i in t0.tolist()])
    set1 = set([tuple(i) for i in t1.tolist()])
    set_inter = set0.intersection(set1)
    set0_idx = {tuple(j):i for i,j in enumerate(t0.tolist())}
    set1_idx = {tuple(j):i for i,j in enumerate(t1.tolist())}
    same_elems_in0 = np.array([set0_idx[i] for i in set_inter])
    same_elems_in1 = np.array([set1_idx[i] for i in set_inter])
    return same_elems_in0, same_elems_in1

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

def advect_nodes_carefully(x0, y0, u0, v0, t, i=0, maxi=100):
    x1 = x0 + u0
    y1 = y0 + v0
    a1 = get_area(x1, y1, t)
    if a1.min() >= 0:
        #print('All areas are positive after advecting nodes carefully', i)
        return x1, y1
    neg_elems = np.nonzero(a1 < 0)[0]
    neg_nodes = np.unique(t[neg_elems].flatten())
    u0[neg_nodes] *= 0.5
    v0[neg_nodes] *= 0.5
    if i >= maxi:
        raise Exception('Max iterations reached in advect_nodes_carefully')
    #print(f'Iteration {i}: {len(neg_elems)} negative elements, {len(neg_nodes)} affected nodes')
    return advect_nodes_carefully(x0, y0, u0, v0, t, i+1, maxi)

def get_mesh_files(idate, lag_dir, mesh_init_file):
    mesh_src_dir = idate.strftime(f'{lag_dir}/%Y')
    mesh_src_file = idate.strftime(f'{mesh_src_dir}/mesh_%Y%m%d.zip')
    if not os.path.exists(mesh_src_file):
        print('Load INIT')
        mesh_src_file = mesh_init_file

    mesh_dst_date = idate + timedelta(1)
    mesh_dst_dir = mesh_dst_date.strftime(f'{lag_dir}/%Y')
    os.makedirs(mesh_dst_dir, exist_ok=True)
    mesh_dst_file = mesh_dst_date.strftime(f'{mesh_dst_dir}/mesh_%Y%m%d.zip')
    return mesh_src_file, mesh_dst_file

def load_data(ifile, mesh_src_file):
    with np.load(ifile) as d:
        u = d['u']
        v = -d['v']
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0
    x0, y0, t0 = MeshFile(mesh_src_file).load(['x', 'y', 't'], as_dict=False)    
    return u, v, Triangulation(x0, y0, t0)

@njit(cache=True, inline="always")
def _cross(ax, ay, bx, by):
    return ax * by - ay * bx


@njit(cache=True)
def _intersection_area(a, b):
    # Translate coordinates to improve numerical precision
    origin_x, origin_y = a[0, 0], a[0, 1]

    p = np.empty((12, 2), dtype=np.float64)
    q = np.empty((12, 2), dtype=np.float64)

    for i in range(3):
        p[i, 0] = a[i, 0] - origin_x
        p[i, 1] = a[i, 1] - origin_y

    c = np.empty((3, 2), dtype=np.float64)
    for i in range(3):
        c[i, 0] = b[i, 0] - origin_x
        c[i, 1] = b[i, 1] - origin_y

    # Scale-dependent tolerance
    scale = 0.0
    for i in range(3):
        for j in range(2):
            scale = max(scale, abs(p[i, j]), abs(c[i, j]))

    if scale == 0.0:
        return 0.0

    eps = 1e-12 * scale * scale

    # Ensure clipping triangle is counterclockwise
    orient = 0.0
    for i in range(3):
        j = (i + 1) % 3
        orient += _cross(
            c[j, 0] - c[i, 0], c[j, 1] - c[i, 1],
            c[(i + 2) % 3, 0] - c[i, 0],
            c[(i + 2) % 3, 1] - c[i, 1]
        )

    # The expression above equals twice the signed triangle area
    if abs(orient) <= eps:
        return 0.0

    if orient < 0.0:
        tmpx, tmpy = c[1, 0], c[1, 1]
        c[1, 0], c[1, 1] = c[2, 0], c[2, 1]
        c[2, 0], c[2, 1] = tmpx, tmpy

    # Reject degenerate subject triangle
    area_a = 0.0
    for i in range(3):
        j = (i + 1) % 3
        area_a += _cross(p[i, 0], p[i, 1],
                         p[j, 0], p[j, 1])

    if abs(area_a) <= eps:
        return 0.0

    n = 3

    # Clip polygon A against each edge of triangle B
    for edge in range(3):
        if n == 0:
            return 0.0

        ex, ey = c[edge, 0], c[edge, 1]
        fx, fy = c[(edge + 1) % 3, 0], c[(edge + 1) % 3, 1]
        dx, dy = fx - ex, fy - ey

        m = 0
        sx, sy = p[n - 1, 0], p[n - 1, 1]
        ds = _cross(dx, dy, sx - ex, sy - ey)
        sin = ds >= -eps

        for i in range(n):
            tx, ty = p[i, 0], p[i, 1]
            dt = _cross(dx, dy, tx - ex, ty - ey)
            tin = dt >= -eps

            # Segment crosses clipping boundary
            if sin != tin:
                denom = ds - dt
                if denom != 0.0:
                    t = ds / denom
                    q[m, 0] = sx + t * (tx - sx)
                    q[m, 1] = sy + t * (ty - sy)
                    m += 1

            if tin:
                q[m, 0] = tx
                q[m, 1] = ty
                m += 1

            sx, sy = tx, ty
            ds, sin = dt, tin

        n = m

        # Swap working buffers
        tmp = p
        p = q
        q = tmp

    # Shoelace formula
    if n < 3:
        return 0.0

    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += _cross(p[i, 0], p[i, 1],
                       p[j, 0], p[j, 1])

    return abs(area) * 0.5


@njit(cache=True)
def triangle_intersection_area(a, b):
    """
    a, b: arrays of shape (3, 2)
    Returns intersection area.
    """
    return _intersection_area(a, b)


@njit(cache=True, parallel=True)
def triangle_intersection_areas(A, B):
    """
    A, B: arrays of shape (N, 3, 2)
    Returns an array of N intersection areas.
    """
    n = A.shape[0]
    result = np.empty(n, dtype=np.float64)

    for i in prange(n):
        result[i] = _intersection_area(A[i], B[i])

    return result


def compute_mapping(tri_a, tri_o, search_dist=50):
    xa, ya, ta = tri_a.x, tri_a.y, tri_a.triangles
    xo, yo, to = tri_o.x, tri_o.y, tri_o.triangles

    same_in_a, same_in_o = get_same_elements(xa[ta], xo[to])
    new_elems = np.ones(tri_o.triangles.shape[0], dtype=bool)
    if same_in_o.size > 0:
        new_elems[same_in_o] = False
    new_elem_idx = np.nonzero(new_elems)[0]

    src2dst = np.column_stack([same_in_a, same_in_o]).tolist()
    weights = [1] * len(src2dst)
    elem_xa = xa[ta].mean(axis=1)
    elem_ya = ya[ta].mean(axis=1)
    elem_xo = xo[to].mean(axis=1)
    elem_yo = yo[to].mean(axis=1)

    tree_a = KDTree(np.column_stack([elem_xa, elem_ya])) 
    area_o = get_area(xo, yo, to)
    for i in new_elem_idx:
        elem_o = np.column_stack((xo[to[i]], yo[to[i]]))
        elem_o_area = area_o[i]
        #print(elem_o)
        near_elems_1 = tree_a.query_ball_point([elem_xo[i], elem_yo[i]], search_dist)
        #print(near_elems_1)
        for j in near_elems_1:
            elem_a = np.column_stack((xa[ta[j]], ya[ta[j]]))
            #print(elem_a)
            inter_area = triangle_intersection_area(elem_o, elem_a)
            #print(inter_area)
            if inter_area > 0:
                src2dst.append((j, i))
                weights.append(inter_area / elem_o_area)
    src2dst = np.array(src2dst)
    weights = np.array(weights)
    return src2dst, weights

def get_area_ratio(tri0, tri_a, tri_o, src2dst, weights):
    area0 = get_area(tri0.x, tri0.y, tri0.triangles)
    area1 = get_area(tri_a.x, tri_a.y, tri_a.triangles)
    area_ratio1 = area1 / area0
    area_ratio7 = np.zeros(tri_o.triangles.shape[0])
    area_weights = np.zeros(tri_o.triangles.shape[0])
    np.add.at(area_ratio7, src2dst[:, 1], area_ratio1[src2dst[:, 0]] * weights)
    np.add.at(area_weights, src2dst[:, 1], weights)
    return area_ratio7 / area_weights

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
