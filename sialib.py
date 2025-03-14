import numpy as np
from matplotlib.tri import Triangulation
import gmsh
from scipy.spatial import ConvexHull

def jacobian(x0, y0, x1, y1, x2, y2):
    return (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)

def get_area(x, y, t):
    return .5*jacobian(x[t][:,0], y[t][:,0], x[t][:,1], y[t][:,1], x[t][:,2], y[t][:,2])

def find_fixed_edges(t, fixed_nodes_idx):
    x_tmp = np.zeros(t.x.size, bool)
    x_tmp[fixed_nodes_idx] = True
    fixed_edges = x_tmp[t.edges[:,0]] & x_tmp[t.edges[:,1]]
    return fixed_edges

def find_movable_edges(t, fixed_nodes_idx):
    return np.nonzero(~find_fixed_edges(t, fixed_nodes_idx))

def remove_short_edges(t, min_edge_length, itern=0, maxiter=1000, fixed_nodes_idx=None):
    if fixed_nodes_idx is not None:
        good_edges = find_movable_edges(t, fixed_nodes_idx)
        t_edges = t.edges[good_edges]
    else:
        t_edges = t.edges

    # Calculate length of each edge
    edge_lengths = np.hypot(np.diff(t.x[t_edges]), np.diff(t.y[t_edges])).flatten()

    # Mask for short edges
    shortest_edge_index = np.argmin(edge_lengths)

    # Stop processing if the shortest edge is longer than the threshold
    if edge_lengths[shortest_edge_index] > min_edge_length or itern >= maxiter:
        return t, itern

    # Get the vertices of this edge
    n1, n2 = t_edges[shortest_edge_index]

    # Calculate midpoint to place the new node
    midpoint_x = (t.x[n1] + t.x[n2]) / 2
    midpoint_y = (t.y[n1] + t.y[n2]) / 2

    new_n = t.x.size
    new_x = np.hstack([t.x, midpoint_x])
    new_y = np.hstack([t.y, midpoint_y])
    new_triangles = np.array(t.triangles)
    new_triangles[new_triangles == n1] = new_n
    new_triangles[new_triangles == n2] = new_n

    edge = sorted(t_edges[shortest_edge_index])
    tri_sorted = np.sort(t.triangles, axis=1)
    remove_triangle_indices = np.nonzero(
        (np.any(tri_sorted == edge[0], axis=1)) &
        (np.any(tri_sorted == edge[1], axis=1))
    )[0]

    good_triangles = np.ones(len(t.triangles), dtype=bool)
    good_triangles[remove_triangle_indices] = False
    new_triangles = new_triangles[good_triangles]

    new_x[n1] = 0 #np.nan
    new_y[n1] = 0 #np.nan
    new_x[n2] = 0 #np.nan
    new_y[n2] = 0 #np.nan
    t = Triangulation(new_x, new_y, new_triangles)
    return remove_short_edges(t, min_edge_length, itern+1, maxiter, fixed_nodes_idx)

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

def prepare_surface_tag(t, model_name="regular_mesh"):
    if not gmsh.is_initialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    # Step 1: Create geometric points
    #valid_node_idx = np.unique(t.triangles)
    valid_node_idx = np.arange(t.x.size)
    point_tags = valid_node_idx + 1
    points = np.column_stack([t.x[valid_node_idx], t.y[valid_node_idx]])

    hull = ConvexHull(points)
    boundary_points = hull.vertices
    for bpi in boundary_points:
        _ = gmsh.model.geo.addPoint(*points[bpi], 0, tag=point_tags[bpi])

    line_tags = []
    for i in range(len(boundary_points)):
        start = point_tags[boundary_points[i]]
        end = point_tags[boundary_points[(i+1)%len(boundary_points)]]
        line_tag = gmsh.model.geo.addLine(start, end)
        line_tags.append(line_tag)

    # Step 3: Create surface
    loop_tag = gmsh.model.geo.addCurveLoop(line_tags)
    surface_tag = gmsh.model.geo.addPlaneSurface([loop_tag])

    # Step 4: Synchronize before meshing
    gmsh.model.geo.synchronize()
    return points, point_tags, surface_tag

def create_mesh_from_points(t, *args, **kwargs):
    points, point_tags, surface_tag = prepare_surface_tag(t, *args, **kwargs)
    node_coords = np.column_stack((points, np.zeros(len(points)))).flatten()
    node_tags = point_tags.tolist() #list(range(1, len(points) + 1))
    gmsh.model.mesh.addNodes(2, 1, node_tags, node_coords)

    # Add triangle elements using the existing triangulation
    element_tags = [[i + 1] for i in range(len(t.triangles))]
    element_node_tags = (t.triangles + 1).tolist()

    # Add elements
    gmsh.model.mesh.addElements(
        2,      # dimension
        surface_tag,  # entity tag (use the created surface tag)
        [2] * len(element_tags),    # element type (2 = triangle)
        element_tags,
        element_node_tags,
    )

def get_triangulation():
    coords = gmsh.model.mesh.getNodes()[1].reshape(-1, 3)[:, :2]
    elem_node_tags = gmsh.model.mesh.getElements(dim=2)[2]
    triangles = elem_node_tags[0].reshape(-1, 3) - 1
    t = Triangulation(coords[:, 0], coords[:, 1], triangles)
    return t

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