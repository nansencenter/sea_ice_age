#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Stitch two nextsim restart files with different meshes.
The first file is remapped to the mesh of the second file.
The output is saved in a npz file in mesh/YYYY/mesh_YYYYMMDD.npz

Usage:
    python stitch_restarts.py <input_file1> <input_file2> <output_dir> [--search-dist 15000] [--cores 1] [--force 0]

Created on Thu Jun 20 10:00:00 2024
"""

import argparse
import os

from matplotlib.tri import Triangulation
from pynextsim import NextsimBin
import numpy as np
from scipy.spatial import cKDTree

from utils import compute_mapping, get_area_ratio


DEFAULT_SEARCH_DIST = 15e3

DEFAULT_CORES = 1


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process two input files and specify an output directory.")
    parser.add_argument('input_file1', type=str, help='Path to the first input file')
    parser.add_argument('input_file2', type=str, help='Path to the second input file')
    parser.add_argument('output_dir', type=str, help='Path to the output directory')
    parser.add_argument('--search-dist', type=float,
            default=DEFAULT_SEARCH_DIST,
            help=f'Search distance for compute_mapping (default: {DEFAULT_SEARCH_DIST})')
    parser.add_argument('--cores', type=int,
            default=DEFAULT_CORES,
            help=f'Number of cores to use (default: {DEFAULT_CORES})')
    parser.add_argument('--force', action="store_true",
            help='Force overwrite of output files')
    return parser.parse_args()


def nextsimbin2tri(restart_file, maskfile='mask.npy', min_x=-2.5e6, max_y=2.1e6, res=20000):
    """ Read Nextsim binary file and mask.
    Return triangulation and node IDs of the masked mesh.

    Parameter:
    ----------
    restart_file: str
        Path to Nextsim binary restart file.
    maskfile: str
        Path to numpy file with mask (2D array of bools).
    min_x, max_y: float
        Coordinates of the upper left corner of the mask.
    res: float
        Resolution of the mask (assumed square pixels).

    Returns:
    ------
    tri: Triangulation
        Triangulation of the masked mesh.
    ids: np.ndarray
        Node IDs of the masked mesh.

    """
    # read raw data
    n = NextsimBin(restart_file)
    tri = Triangulation(n.mesh_info.nodes_x, n.mesh_info.nodes_y, n.mesh_info.indices)
    ids = n.mesh_info.get_var('id')

    # read mask
    mask = np.load(maskfile)
    # convert element coordinates to mask indices
    x_el = tri.x[tri.triangles].mean(axis=1)
    y_el = tri.y[tri.triangles].mean(axis=1)
    cols_el = np.clip((x_el - min_x) / res, 0, mask.shape[1]-1).astype(int)
    rows_el = np.clip((max_y - y_el) / res, 0, mask.shape[0]-1).astype(int)
    # read mask values for each element
    el_mask = mask[rows_el, cols_el]
    # subset triangulation by masked elements
    sub_tri = tri.triangles[el_mask]
    # regenerate triangulation by excluding unused nodes
    uniq_nodes, uniq_inv = np.unique(sub_tri, return_inverse=true)
    newx = tri.x[uniq_nodes]
    newy = tri.y[uniq_nodes]
    newt = uniq_inv.reshape(-1,3)
    new_ids = ids[uniq_nodes]
    # create new triangulation
    tri = Triangulation(newx, newy, newt)
    ids = new_ids
    return tri, ids


def get_tri_a_from_nextsim(tri_0, ids_0, tri_o, ids_o, k=10):
    """ Create triangulation tri_a with the same connectivity as tri_0,
    but with node coordinates from tri_o where possible.
    For nodes in tri_0 that are not present in tri_o, find the nearest
    neighbor in tri_0 that is present in tri_o, and use the same offset
    from that neighbor to assign coordinates.

    Parameters:
    -----------
    tri_0: Triangulation
        Triangulation of the first mesh (to be remapped).
    ids_0: np.ndarray
        Node IDs of the first mesh.
    tri_o: Triangulation
        Triangulation of the second mesh (target mesh).
    ids_o: np.ndarray
        Node IDs of the second mesh.
    k: int
        Number of nearest neighbors to search for each lost node.
    
    Returns:
    --------
    tri_a: Triangulation
        Triangulation with the same connectivity as tri_0, but with
        node coordinates from tri_o where possible.

    """
    # find common nodes and create new nodes with coordinates from the target triangulation
    new_x = np.zeros_like(tri_0.x) + np.nan
    new_y = np.zeros_like(tri_0.y) + np.nan
    _, idx0, idx1 = np.intersect1d(ids_0, ids_o, return_indices=True)
    new_x[idx0] = tri_o.x[idx1]
    new_y[idx0] = tri_o.y[idx1]

    # for lost nodes, find nearest neighbor in tri_0 that is present in tri_o
    # and use the same offset from that neighbor to assign coordinates
    points = np.column_stack((tri_0.x, tri_0.y))
    tree = cKDTree(points)
    lost_node_ids = list(set(ids_0) - set(ids_o))
    for lost_node_id in lost_node_ids:
        lost_node_index = np.nonzero(ids_0 == lost_node_id)[0]
        lost_node_x = tri_0.x[lost_node_index]
        lost_node_y = tri_0.y[lost_node_index]

        # increase k until we find at least one common neighbor
        for k_factor in [1, 2, 5, 10]:
            # find k nearest neighbors in tri_0
            distances, neib_indices = tree.query(points[lost_node_index], k=k * k_factor)
            neib_ids = ids_0[neib_indices]
            # find which of these neighbors are present in tri_o
            common_neib_ids, idx0, idx1 = np.intersect1d(neib_ids, ids_o, assume_unique=True, return_indices=True)
            # if common neighbors found, break
            if len(common_neib_ids) > 0:
                break

        # find the nearest neighbor among the common neighbors
        nearest_neib_subindex = np.argmin(distances[:, idx0])
        # get the index of this neighbor in tri_0 and tri_o
        nearest_neib_index0 = neib_indices[:, idx0[nearest_neib_subindex]]
        nearest_neib_index1 = idx1[nearest_neib_subindex]
        # get coordinates of this neighbor in tri_0
        nn_node_x0 = tri_0.x[nearest_neib_index0]
        nn_node_y0 = tri_0.y[nearest_neib_index0]
        # compute offset from neighbor to lost node
        dx = lost_node_x - nn_node_x0
        dy = lost_node_y - nn_node_y0
        # assign new coordinates to lost node using the same offset from the neighbor in tri_o
        new_x[lost_node_index] = tri_o.x[nearest_neib_index1] + dx
        new_y[lost_node_index] = tri_o.y[nearest_neib_index1] + dy

    tri_a = Triangulation(new_x, new_y, triangles=tri_0.triangles)
    return tri_a


def get_dst_file(input_file2, output_dir):
    """ Create output file path based on the date in input_file2."""
    date_str = os.path.basename(input_file2).split('.')[0].split('_')[1].split('T')[0]
    yyyy_str = date_str[:4]
    mesh_dst_dir = f'{output_dir}/{yyyy_str}'
    os.makedirs(mesh_dst_dir, exist_ok=True)
    mesh_dst_file = f'{mesh_dst_dir}/mesh_{date_str}.npz'
    return mesh_dst_file


def stitch_restart_pair(input_file1, input_file2, output_dir,
        force=False, search_dist=DEFAULT_SEARCH_DIST,
        cores=DEFAULT_CORES):

    print(f"Input File 1: {input_file1}")
    print(f"Input File 2: {input_file2}")
    print(f"Output Directory: {output_dir}")
    print(f"Search Distance: {search_dist}")
    print(f"Cores: {cores}")
    print(f"Force Overwrite: {force}")
    mesh_dst_file = get_dst_file(input_file2, output_dir)
    if os.path.exists(mesh_dst_file) and not force:
        print(f"Output file {mesh_dst_file} already exists. Use --force 1 to overwrite.")
        return
    
    tri_0, ids_0 = nextsimbin2tri(input_file1)
    tri_o, ids_o = nextsimbin2tri(input_file2)
    tri_a = get_tri_a_from_nextsim(tri_0, ids_0, tri_o, ids_o)
    src2dst, weights = compute_mapping(tri_a, tri_o, search_dist, cores=cores)
    area_ratio = get_area_ratio(tri_0, tri_a, tri_o, src2dst, weights)
    print(f"Saving {mesh_dst_file}")
    np.savez(mesh_dst_file, x=tri_o.x, y=tri_o.y, t=tri_o.triangles, src2dst=src2dst, weights=weights, ar=area_ratio, ids=ids_o)


def main():
    args = parse_arguments()
    stitch_restart_pair(**vars(args))


if __name__ == "__main__":
    main()
