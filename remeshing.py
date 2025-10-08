import numpy as np
from matplotlib.tri import Triangulation

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

def find_fixed_edges(t, fixed_nodes_idx, strict=False):
    x_tmp = np.zeros(t.x.size, bool)
    x_tmp[fixed_nodes_idx] = True
    if strict:
        fixed_edges = x_tmp[t.edges[:,0]] | x_tmp[t.edges[:,1]]
    else:
        fixed_edges = x_tmp[t.edges[:,0]] & x_tmp[t.edges[:,1]]
    return fixed_edges

def find_movable_edges(t, fixed_nodes_idx, strict=False):
    return np.nonzero(~find_fixed_edges(t, fixed_nodes_idx, strict))

def collapse_short_edges(t, min_edge_length, itern=0, maxiter=1000, fixed_nodes_idx=None):
    if fixed_nodes_idx is not None:
        good_edges = find_movable_edges(t, fixed_nodes_idx, strict=True)
        t_edges = t.edges[good_edges]
    else:
        t_edges = t.edges

    edge_lengths = np.hypot(np.diff(t.x[t_edges]), np.diff(t.y[t_edges])).flatten()
    shortest_edge_index = np.argmin(edge_lengths)
    if edge_lengths[shortest_edge_index] > min_edge_length or itern >= maxiter:
        return t, itern
    t = collapse_edge(t, *t_edges[shortest_edge_index])
    return collapse_short_edges(t, min_edge_length, itern+1, maxiter, fixed_nodes_idx)

def collapse_edge(t, n1, n2):
    midpoint_x = (t.x[n1] + t.x[n2]) / 2
    midpoint_y = (t.y[n1] + t.y[n2]) / 2

    new_n = t.x.size
    new_x = np.hstack([t.x, midpoint_x])
    new_y = np.hstack([t.y, midpoint_y])
    new_triangles = np.array(t.triangles)
    new_triangles[new_triangles == n1] = new_n
    new_triangles[new_triangles == n2] = new_n

    good_triangles = np.all(t.triangles != n1, axis=1) | np.all(t.triangles != n2, axis=1)
    new_triangles = new_triangles[good_triangles]

    new_x[n1] = 0
    new_y[n1] = 0
    new_x[n2] = 0
    new_y[n2] = 0
    return Triangulation(new_x, new_y, new_triangles)

def get_all_angles(t, x, y):
    angles = np.zeros(t.shape)
    for i in range(3):
        vec1list = np.column_stack([x[t[:, i]] - x[t[:, (i+2)%3]], y[t[:,i]] - y[t[:, (i+2)%3]]])
        vec2list = np.column_stack([x[t[:, (i+1)%3]] - x[t[:, (i+2)%3]], y[t[:, (i+1)%3]] - y[t[:, (i+2)%3]]])
        dot_product = np.sum(vec1list * vec2list, axis=1)
        norm1 = np.linalg.norm(vec1list, axis=1)
        norm2 = np.linalg.norm(vec2list, axis=1)
        cos_angle = dot_product / (norm1 * norm2)
        angles[:, (i+2)%3] = np.degrees(np.arccos(cos_angle))
    return angles

def get_min_angle_edge(x, y, t, min_angle=10, fixed_nodes_idx=None):
    angles = get_all_angles(t, x, y)
    tri_with_smal_angles = np.nonzero(np.any(angles < min_angle, axis=1))[0]
    if tri_with_smal_angles.size == 0:
        return None
    min_angle_corner = np.argmin(angles[tri_with_smal_angles], axis=1)
    min_angle_edge1 = (min_angle_corner + 1)%3
    min_angle_edge2 = (min_angle_corner + 2)%3
    min_angle_edges = np.column_stack([
        t[tri_with_smal_angles, min_angle_edge1],
        t[tri_with_smal_angles, min_angle_edge2]
    ])
    if fixed_nodes_idx is None:
        return min_angle_edges[0]
    x_tmp = np.zeros(x.size, bool)
    x_tmp[fixed_nodes_idx] = True
    fixed_edges = x_tmp[min_angle_edges[:,0]] | x_tmp[min_angle_edges[:,1]]
    min_angle_edges = min_angle_edges[~fixed_edges]
    if min_angle_edges.size  == 0:
        return None
    return min_angle_edges[0]

def collapse_small_angle_edges(t, edge_angle_threshold=5, niter=0, max_iter=100, fixed_nodes_idx=None):
    min_angle_edge = get_min_angle_edge(t.x, t.y, t.triangles, edge_angle_threshold, fixed_nodes_idx)
    if min_angle_edge is None or niter >= max_iter:
        return t, niter
    t = collapse_edge(t, *min_angle_edge)
    return collapse_small_angle_edges(t, edge_angle_threshold, niter+1, max_iter, fixed_nodes_idx=fixed_nodes_idx)

def remove_long_edges(t, max_edge_length, itern=0, maxiter=1000, fixed_nodes_idx=None):
    if fixed_nodes_idx is not None:
        good_edges = find_movable_edges(t, fixed_nodes_idx)
        t_edges = t.edges[good_edges]
    else:
        t_edges = t.edges

    edge_lengths = np.hypot(np.diff(t.x[t_edges]), np.diff(t.y[t_edges])).flatten()
    longest_edge_index = np.argmax(edge_lengths)

    if edge_lengths[longest_edge_index] < max_edge_length or itern >= maxiter:
        return t, itern

    n1, n2 = t_edges[longest_edge_index]
    midpoint_x = (t.x[n1] + t.x[n2]) / 2
    midpoint_y = (t.y[n1] + t.y[n2]) / 2

    new_n = t.x.size
    new_x = np.hstack([t.x, midpoint_x])
    new_y = np.hstack([t.y, midpoint_y])
    new_triangles = np.array(t.triangles)

    edge = sorted(t_edges[longest_edge_index])
    tri_sorted = np.sort(t.triangles, axis=1)
    remove_triangle_indices = sorted(np.nonzero(
        (np.any(tri_sorted == edge[0], axis=1)) &
        (np.any(tri_sorted == edge[1], axis=1))
    )[0], reverse=True)

    add_triangles = []
    for ti in remove_triangle_indices:
        for e in t_edges[longest_edge_index]:
            add_triangle = np.array(t.triangles[ti])
            add_triangle[add_triangle == e] = new_n
            add_triangles.append(add_triangle)
    add_triangles = np.array(add_triangles)
    new_triangles = list(new_triangles)
    for idx in remove_triangle_indices:
        new_triangles.pop(idx)
    new_triangles = np.vstack([new_triangles, add_triangles])
    t = Triangulation(new_x, new_y, new_triangles)
    return remove_long_edges(t, max_edge_length, itern+1, maxiter, fixed_nodes_idx)

def find_triangle(x, y, t, point):
    # Reshape to get triangle vertices as (n_triangles, 3, 2)
    vertices = np.dstack([x[t], y[t]])

    # Calculate vectors
    v0 = vertices[:, 1] - vertices[:, 0]
    v1 = vertices[:, 2] - vertices[:, 0]

    # Broadcast the point to calculate v2 for all triangles at once
    point_array = np.array(point)
    v2 = point_array - vertices[:, 0]

    # Compute dot products
    d00 = np.sum(v0 * v0, axis=1)
    d01 = np.sum(v0 * v1, axis=1)
    d11 = np.sum(v1 * v1, axis=1)
    d20 = np.sum(v2 * v0, axis=1)
    d21 = np.sum(v2 * v1, axis=1)

    # Calculate denominators
    denom = d00 * d11 - d01 * d01

    # Avoid division by zero
    valid_denom = denom != 0

    # Initialize barycentric coordinates
    v = np.full_like(denom, np.inf)
    w = np.full_like(denom, np.inf)

    # Calculate only where denominator is valid
    v[valid_denom] = (d11[valid_denom] * d20[valid_denom] - d01[valid_denom] * d21[valid_denom]) / denom[valid_denom]
    w[valid_denom] = (d00[valid_denom] * d21[valid_denom] - d01[valid_denom] * d20[valid_denom]) / denom[valid_denom]
    u = 1.0 - v - w

    # Point is inside triangle if all barycentric coordinates are positive
    inside = (u > 0) & (v > 0) & (w > 0)

    # Return the index of the first triangle containing the point, or -1 if none
    containing_triangles = np.where(inside)[0]
    return containing_triangles[0] if len(containing_triangles) > 0 else -1

def find_triangles_for_points(x, y, t, points):
    return np.array([find_triangle(x, y, t, point) for point in points])

def order_cclockwise(x, y):
    """
    Order the points defined by x and y coordinates in counter clockwise order.

    Parameters:
    -----------
    x : array-like
        x-coordinates of the points
    y : array-like
        y-coordinates of the points

    Returns:
    --------
    x_ordered : array-like
        x-coordinates of the points in clockwise order
    y_ordered : array-like
        y-coordinates of the points in clockwise order
    """
    # Find centroid
    centroid_x = np.mean(x)
    centroid_y = np.mean(y)

    # Calculate angles from centroid to each point
    angles = np.arctan2(y - centroid_y, x - centroid_x)

    # Sort points by angle (clockwise)
    sorted_indices = np.argsort(angles)

    # Return sorted coordinates
    return sorted_indices

def remove_filpped_elements(x, y, t, niter=0, max_iter=1000):
    area = get_area(x, y, t)
    min_area_i = np.argmin(area)
    if area[min_area_i] > 0 or niter >= max_iter:
        return x, y, t, niter

    neg_tri_x = x[t[min_area_i]].flatten()
    neg_tri_y = y[t[min_area_i]].flatten()

    potential_trouble_points = np.column_stack([neg_tri_x, neg_tri_y])
    tri_i = find_triangles_for_points(x, y, t, potential_trouble_points)

    if np.all(tri_i == -1):
        reorder_index = Triangulation(x[t[min_area_i]], y[t[min_area_i]]).triangles[0]
        t[min_area_i] = t[min_area_i][reorder_index]
        return x, y, t, niter

    trouble_node = t[min_area_i][tri_i > 0][0]
    edges = np.unique(np.sort(np.vstack([t[:, 0:2], t[:, 1:3], np.column_stack([t[:, 2], t[:, 0]])]), axis=1)[:, ::-1], axis=0)
    trouble_edges = np.nonzero(np.any(edges == trouble_node, axis=1))
    neib_nodes = np.unique(edges[trouble_edges])
    new_node_x = x[neib_nodes].mean()
    new_node_y = y[neib_nodes].mean()

    bad_triangs = np.any(t == trouble_node, axis=1)
    new_t = np.array(t)
    new_x = np.hstack([x, new_node_x])
    new_y = np.hstack([y, new_node_y])
    add_triangles = np.array(t[bad_triangs])
    add_triangles[add_triangles == trouble_node] = new_x.size - 1
    new_t[bad_triangs] = add_triangles
    return remove_filpped_elements(new_x, new_y, new_t, niter + 1, max_iter)

def find_hanging_elements(x, y, t):
    edges = np.unique(np.sort(np.vstack([t[:, 0:2], t[:, 1:3], np.column_stack([t[:, 2], t[:, 0]])]), axis=1)[:, ::-1], axis=0)
    edge_count = np.zeros(x.size)
    np.add.at(edge_count, edges, 1)
    suspicious_node_ids = np.nonzero(edge_count==2)[0]
    potential_trouble_points = np.column_stack([x[suspicious_node_ids], y[suspicious_node_ids]])
    tri_i = find_triangles_for_points(x, y, t, potential_trouble_points)
    bad_nodes = suspicious_node_ids[tri_i >= 0]
    bad_elements = np.array([np.nonzero((t == bad_node).any(axis=1))[0] for bad_node in bad_nodes]).flatten()
    return bad_elements

def remove_hanging_elements(tri):
    bad_elements = find_hanging_elements(tri.x, tri.y, tri.triangles)
    if bad_elements.size == 0:
        return tri
    good_tri_mask = np.ones(tri.triangles.shape[0], dtype=bool)
    good_tri_mask[bad_elements] = False
    t_out = np.array(tri.triangles)
    t_out = t_out[good_tri_mask]
    tri_out = Triangulation(tri.x, tri.y, t_out)
    return tri_out

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

def laplace_smooth(tri, iterations=1, factor=0.5, fixed_nodes=None):
    """
    Apply Laplace smoothing to a triangulation mesh.

    Parameters:
    -----------
    tri : Triangulation
        The triangulation to smooth
    iterations : int
        Number of smoothing iterations
    factor : float
        Relaxation factor (0-1), how much to move toward average position
    fixed_nodes : array-like
        Indices of nodes that should not be moved

    Returns:
    --------
    new_tri : Triangulation
        Smoothed triangulation
    """
    # Make a copy of the triangulation to work with
    x = tri.x.copy()
    y = tri.y.copy()
    triangles = np.array(tri.triangles)

    # Create a mask for fixed nodes if provided
    mask = np.ones(len(x), dtype=bool)
    if fixed_nodes is not None:
        mask[fixed_nodes] = False

    # Get the adjacency information
    edges = {}
    for t in triangles:
        for i, j in [(0, 1), (1, 2), (2, 0)]:
            a, b = t[i], t[j]
            if a not in edges:
                edges[a] = []
            if b not in edges:
                edges[b] = []
            if b not in edges[a]:
                edges[a].append(b)
            if a not in edges[b]:
                edges[b].append(a)

    # Perform smoothing iterations
    for _ in range(iterations):
        new_x = x.copy()
        new_y = y.copy()

        # For each node, compute new position
        for i in range(len(x)):
            if mask[i] and i in edges:  # Skip fixed nodes
                neighbors = edges[i]
                if neighbors:
                    # Compute average position of neighbors
                    avg_x = np.mean([x[j] for j in neighbors])
                    avg_y = np.mean([y[j] for j in neighbors])

                    # Move toward average position by factor
                    new_x[i] = x[i] + factor * (avg_x - x[i])
                    new_y[i] = y[i] + factor * (avg_y - y[i])

        # Update positions
        x, y = new_x, new_y

    # Create new triangulation with smoothed points
    new_tri = Triangulation(x, y, triangles)
    return new_tri

def sort_all_nodes_ccw(x, y, t):
    centroid_x = np.mean(x[t], axis=1)[None].T
    centroid_y = np.mean(y[t], axis=1)[None].T
    angles = np.arctan2(y[t] - centroid_y, x[t] - centroid_x)
    sorted_indices = np.argsort(angles, axis=1)
    rows = np.column_stack([np.arange(t.shape[0])]*3)
    cols = sorted_indices
    t7rot = t[rows, cols]
    return t7rot

def remove_element(x, y, t, fixed_nodes_idx, neg_tri_i):
    neg_tri = t[neg_tri_i]
    neg_tri_i_edges = np.array([[neg_tri[i], neg_tri[(i+1)%3].item()] for i in range(3)])
    edge_lengths = np.hypot(x[neg_tri_i_edges[:,0]] - x[neg_tri_i_edges[:,1]], y[neg_tri_i_edges[:,0]] - y[neg_tri_i_edges[:,1]])
    valid_edges = np.nonzero([(edge[0] not in fixed_nodes_idx) or (edge[1] not in fixed_nodes_idx) for edge in neg_tri_i_edges])[0]
    min_valid_edge = neg_tri_i_edges[valid_edges[np.argmin(edge_lengths[valid_edges])]]
    if min_valid_edge[0] in fixed_nodes_idx:
        good_node = min_valid_edge[0].item()
        bad_node = min_valid_edge[1].item()
    else:
        good_node = min_valid_edge[1].item()
        bad_node = min_valid_edge[0].item()

    new_triangles = np.array(t)
    new_triangles[new_triangles == bad_node] = good_node
    new_triangles = new_triangles[np.sum(new_triangles == good_node, axis=1) < 2]
    new_x = np.array(x)
    new_y = np.array(y)
    new_x[bad_node] = 0
    new_y[bad_node] = 0
    return new_x, new_y, new_triangles

def replace_negative_element_with_two(x, y, t, neg_area_i):
    neg_tri_x = x[t[neg_area_i]].flatten()
    neg_tri_y = y[t[neg_area_i]].flatten()
    potential_trouble_points = np.column_stack([neg_tri_x, neg_tri_y])
    tri_i = find_triangles_for_points(x, y, t, potential_trouble_points)
    contain_tri_i = tri_i[tri_i > 0][0]

    trouble_node = t[neg_area_i][tri_i > 0][0]
    flipped_tri = t[neg_area_i]
    contain_tri = t[contain_tri_i]
    common_edge = np.intersect1d(flipped_tri, contain_tri)
    opposite_node = list(set(list(contain_tri)) - set(list(flipped_tri)))[0]
    new_tri1 = np.array([common_edge[0], opposite_node, trouble_node])
    new_tri2 = np.array([common_edge[1], opposite_node, trouble_node])
    new_tri1 = new_tri1[order_cclockwise(x[new_tri1], y[new_tri1])]
    new_tri2 = new_tri2[order_cclockwise(x[new_tri2], y[new_tri2])]
    new_triangles = np.array(t)
    good_new_tri = np.ones(new_triangles.shape[0], dtype=bool)
    good_new_tri[neg_area_i] = False
    good_new_tri[contain_tri_i] = False
    new_triangles = new_triangles[good_new_tri]
    new_triangles = np.vstack([new_triangles, new_tri1, new_tri2])
    return new_triangles

def get_min_area_for_valid_elemenst(x, y, t, area, fixed_nodes_idx=None):
    if fixed_nodes_idx is None:
        min_area_i = np.argmin(area)
        return area[min_area_i], min_area_i
    fixed_x_flag = np.zeros(x.size, dtype=bool)
    fixed_x_flag[fixed_nodes_idx] = True
    fixed_elem_flag = fixed_x_flag[t].sum(axis=1) == 3
    area_min = area[~fixed_elem_flag].min()
    min_area_i = np.nonzero(area == area_min)[0][0]
    return area_min, min_area_i

def remove_small_elements(tri, fixed_nodes_idx, min_area=20, max_iter=100):
    x, y, t = tri.x, tri.y, tri.triangles
    area = get_area(x, y, t)
    area_min, min_area_i = get_min_area_for_valid_elemenst(x, y, t, area, fixed_nodes_idx)

    n = 0
    while area_min < min_area:
        #print(n, area_min)
        x, y, t = remove_element(x, y, t, fixed_nodes_idx, min_area_i)
        area = get_area(x, y, t)

        if area.min() < 0 and area.min() < area_min:
            #print('Replace with two')
            t = replace_negative_element_with_two(x, y, t, np.argmin(area))
            area = get_area(x, y, t)
        area_min, min_area_i = get_min_area_for_valid_elemenst(x, y, t, area, fixed_nodes_idx)

        n += 1
        if n >= max_iter:
            break
    return Triangulation(x, y, t), n