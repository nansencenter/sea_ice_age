import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RectBivariateSpline
from tqdm import tqdm

from lmsiage.utils import fill_gaps
from lmsiage.mesh_file import MeshFile
from lmsiage.zarr_index_tools import update_index_for_dir, cleanup_missing_files
from lmsiage.zarr_index import get_file_arrays


sid_dir = 'OSISAF_ice_drift_CDR_postproc'
sia_dir = 'NERSC_arctic25km_sea_ice_age_v2p1/zarr'
mesh_dir = f'{sia_dir}/mesh'
age_dir = f'{sia_dir}/age'

update_index_for_dir(mesh_dir)
cleanup_missing_files()

# Load all inital values
mesh_init_file = 'mesh_arctic_ease_25km_max7.npz'
xc = np.load(mesh_init_file)['xc']
yc = np.load(mesh_init_file)['yc']
landmask = np.load(mesh_init_file)['landmask']

min_lm = 0.5
force = False

def interpolate(mesh_file):
    basename = os.path.basename(mesh_file)
    mesh_date = basename.split('.')[0].split('_')[-1]
    mesh_year = mesh_date[:4]

    mf = MeshFile(mesh_file)
    mesh_arrays = get_file_arrays(mesh_file)
    if not force and 'sic' in mesh_arrays:
        return

    x, y, t = mf.load(['x', 'y', 't'], as_dict=False)
    file_mask = f'{sid_dir}/{mesh_year}/ice_drift_nh_ease*{mesh_date}1200.nc.npz'
    try:
        sic_src_file = glob.glob(file_mask)[0]
    except:
        raise ValueError(f'Cannot find {file_mask}')
    try:
        with np.load(sic_src_file) as data:
            cgrd = data['c']
    except:
        print(f'Cannot load c from {sic_src_file}')
        raise ValueError

    # NB: Major change from 2.1
    # Extrapolation of SIC into coast
    # cgrd[np.isnan(cgrd)] = 0
    cgrd_f = fill_gaps(cgrd, np.isnan(cgrd), 10)
    cgrd_f[np.isnan(cgrd_f)] = 0

    try:
        c = RectBivariateSpline(xc, yc, cgrd_f[::-1], kx=1, ky=1)(y[t].mean(axis=1), x[t].mean(axis=1), grid=False)
    except:
        raise ValueError(f'Fail to interpolate SIC in {mesh_file}')
    try:
        lm = RectBivariateSpline(xc, yc, landmask[::-1], kx=1, ky=1)(y[t].mean(axis=1), x[t].mean(axis=1), grid=False)
    except:
        raise ValueError(f'Fail to interpolate LM in {mesh_file}')
    
    c[lm > min_lm] = 0

    data = {'sic': c}
    if 'sic' in mf.read_names():
        mode = 'o'
    else:
        mode = 'a'
    mf.save(data, mode=mode)

mesh_files = sorted(glob.glob(f'{mesh_dir}/202[5,6]/mesh_*zip'))
print(len(mesh_files), mesh_files[0], mesh_files[-1])

for mesh_file in tqdm(mesh_files):
    interpolate(mesh_file)
