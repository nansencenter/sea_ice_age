import glob
import os
from datetime import datetime

from bamg_remesher import BamgRemesher, BamgOpts
from matplotlib.tri import Triangulation
import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter
from tqdm import tqdm
import yaml

from lmsiage.mesh_file import MeshFile


def load_config(config_file='move_mesh.yaml'):
    """Load configuration from YAML file."""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config

# Usage
#if __name__ == '__main__':
config = load_config()
print(config)


# Load mesh initialization data
with np.load(config['mesh_init_file']) as data:
    xc = data['xc']
    yc = data['yc']
    mask = data['mask']
    x0 = data['x']
    y0 = data['y']
    t0 = data['t']
mask_f = gaussian_filter(mask.astype(float), 1, truncate=1)


# Create a grid of points and find which ones are inside the triangulation
grid_resolution = int(config['grid_resolution'])
xgrid, ygrid = np.meshgrid(
    np.arange(xc.min() - grid_resolution, xc.max() - grid_resolution, grid_resolution),
    np.arange(yc.min() - grid_resolution, yc.max() - grid_resolution, grid_resolution),
)
tri0 = Triangulation(x0, y0, t0)
tf0 = tri0.get_trifinder()
idx0 = tf0(xgrid, ygrid)
mask = maximum_filter(idx0 >= 0, size=20 // grid_resolution)
xgrid = xgrid[mask]
ygrid = ygrid[mask]


# Generate start mesh (if it doesn't exist) and save it
lag_dir = config['lag_dir']
mesh_init_file = f'{lag_dir}/1991/mesh_19910101.zip'
if not os.path.exists(mesh_init_file):
    opts = BamgOpts(hmin=15, hmax=25)
    remesher = BamgRemesher()
    result = remesher.remesh(t0, x0, y0, opts)
    t1 = result['triangles']
    x1 = result['x_coords']
    y1 = result['y_coords']
    MeshFile(mesh_init_file).save({'x': x1, 'y': y1, 't': t1})


# Find forcing files and their dates
sid_dir = config['sid_dir']
sid_mask = config['sid_mask']
ifiles = sorted(glob.glob(f'{sid_dir}/{sid_mask}'))
idates = [datetime.strptime(os.path.basename(ifile).split('-')[-1].split('.')[0], '%Y%m%d%H%M%S')
          for ifile in ifiles]
print(len(ifiles), idates[0], idates[-1])


# MOVE MESH
for ifile, idate in tqdm(zip(ifiles, idates), total=len(ifiles)):
    mesh_src_file, mesh_dst_file = get_mesh_files(idate, lag_dir, mesh_init_file)    
    #print(f'Processing {idate} with source mesh {mesh_src_file} and destination mesh {mesh_dst_file}')
    if os.path.exists(mesh_dst_file) and not force:
        continue
    
    # Load U, V from forcing
    with np.load(ifile) as d:
        u = d['u']
        v = -d['v']
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0

    # Load previous mesh
    x0, y0, t0 = MeshFile(mesh_src_file).load(['x', 'y', 't'], as_dict=False)
    tri0 = Triangulation(x0, y0, t0)


    u0, v0, fixed_nodes_idx = resample_velocities(u, v, mask_f, xc, yc, tri0.x, tri0.y, max_mask_fix_node, max_mask_zero_speed)
    x_a, y_a = advect_nodes_carefully(tri0.x, tri0.y, u0, v0, tri0.triangles)
    tri_a = Triangulation(x_a, y_a, tri0.triangles)

    result = remesher.remesh(tri_a.triangles, tri_a.x, tri_a.y, opts)

    t_o = result['triangles']
    x_o = result['x_coords']
    y_o = result['y_coords']    
    tri_o = Triangulation(x_o, y_o, t_o)
    src2dst, weights = compute_mapping_grid(tri_a, tri_o)
    area_ratio = get_area_ratio(tri0, tri_a, tri_o, src2dst, weights)

    #print(f'Save mesh file to {mesh_dst_file}')
    MeshFile(mesh_dst_file).save(dict(
        x=tri_o.x,
        y=tri_o.y,
        t=tri_o.triangles,
        src2dst=src2dst,
        weights=weights,
        ar=area_ratio))
