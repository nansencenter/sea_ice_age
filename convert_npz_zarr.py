from datetime import datetime
import os
import glob
import numpy as np
from lmsiage.mesh_file import MeshFile
from tqdm import tqdm

idir = 'NERSC_arctic25km_sea_ice_age_v2p1'
odir = 'NERSC_arctic25km_sea_ice_age_v2p1/zarr/mesh'
force = False

# load and save mesh files
for year in tqdm(range(2000, 2025)):
    mesh_files = sorted(glob.glob(os.path.join(f'{idir}/mesh/{year}/*.npz')))
    for mesh_file in tqdm(mesh_files):
        mesh_file_date = mesh_file.split('_')[-1].split('.npz')[0]
        mesh_file_date = datetime.strptime(mesh_file_date, '%Y%m%d')
        # output file
        odir_year = f'{odir}/{mesh_file_date.year}'
        ofile = f'{odir_year}/mesh_{mesh_file_date.strftime("%Y%m%d")}.zip'
        if not force and os.path.exists(ofile) and 'x' in MeshFile(ofile).read_names():
            continue

        # find all source files
        try:
            with np.load(mesh_file, allow_pickle=True) as f:
                mesh_data = dict(**f)
        except Exception as e:
            print(f'Error loading mesh file: {mesh_file}, {e}')
            continue

        # save to zarr
        os.makedirs(odir_year, exist_ok=True)
        mf = MeshFile(ofile)
        mf.save(mesh_data)


