import os
import glob

from tqdm import tqdm
import numpy as np
from pynextsim import NextsimBin
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.tri import Triangulation
import zarr

from lmsiage.mesh_file import MeshFile

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
    sic = n.get_var('M_conc') * 100 # for compatibility with previous code...
    sit = n.get_var('M_thick')

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
    uniq_nodes, uniq_inv = np.unique(sub_tri.flatten(), return_inverse=True)
    newx = tri.x[uniq_nodes]
    newy = tri.y[uniq_nodes]
    newt = uniq_inv.reshape(-1,3)
    new_ids = ids[uniq_nodes]
    new_sic = sic[el_mask]
    new_sit = sit[el_mask]
    # create new triangulation
    tri = Triangulation(newx, newy, newt)
    return tri, new_ids, new_sic, new_sit


res_dir = 'restarts'
mesh_npz_dir = 'mesh_npz'
mesh_dir = 'mesh'
years = [1992]


for year in years:
    mesh_npz_files = sorted(glob.glob(f'{mesh_npz_dir}/{year}/mesh*.npz'))
    if year == 1991:
        mesh_npz_files = mesh_npz_files[240:]
    for mesh_npz_file in tqdm(mesh_npz_files):
        mesh_date = os.path.basename(mesh_npz_file).split('_')[1].split('.')[0]
        mesh_yyyy = mesh_date[:4]
        restart_date = (pd.Timestamp(mesh_date) + pd.Timedelta(days=0)).strftime('%Y%m%d')
        restart_file = f'{res_dir}/{restart_date}/inputs/field_{restart_date}T000000Z.bin'
        dst_dir = f'{mesh_dir}/{mesh_yyyy}'
        dst_file = f'{dst_dir}/{mesh_date}.zip'
        if os.path.exists(dst_file):
            continue
        mesh_file = MeshFile(dst_file)
        os.makedirs(dst_dir, exist_ok = True)
        tri, ids, sic, sit = nextsimbin2tri(restart_file)

        with np.load(mesh_npz_file) as f:
            x = f['x']
            y = f['y']
            ids = f['ids']

            t = f['t']
            ar = f['ar']

            src2dst = f['src2dst']
            weights = f['weights']
        
        data = {
            'x': x.astype(np.float32),
            'y': y.astype(np.float32),
            'i': ids,
            't': t,
            'ar': ar.astype(np.float32),
            'src2dst': src2dst,
            'w': weights.astype(np.float32),
            'sic': sic.astype(np.float32),
            'sit': sit.astype(np.float32),
        }
        mesh_file.save(data)
        print(f'Saved {dst_file}')