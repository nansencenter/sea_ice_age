import os

from tqdm import tqdm

from lmsiage.uncertainty import ComputeAgeUncertainty
from lmsiage.zarr_index_tools import cleanup_missing_files, files_with_arrays, update_index_for_dir, files_missing_array

# set up dir names
force = False
sia_dir = 'NERSC_arctic25km_sea_ice_age_v2p1/zarr'
mesh_dir = f'{sia_dir}/mesh'
age_dir = f'{sia_dir}/age'
unc_dir = f'{sia_dir}/uncert'

# reindex files
update_index_for_dir(unc_dir)
cleanup_missing_files()

# find unprocessed uncertainty files
ready_unc_files = files_with_arrays(['unc_sic', 'unc_sid'])
ready_unc_files = [f for f in ready_unc_files if os.path.basename(f).split('_')[1].split('.')[0] > '19910914']
if force:
    unproc_unc_files = ready_unc_files
else:
    unproc_unc_files = files_missing_array('unc_age')
    unproc_unc_files = sorted(list(set(ready_unc_files).intersection(set(unproc_unc_files))))
print(unproc_unc_files[0], unproc_unc_files[-1], len(unproc_unc_files))    

# compute age uncertainty for unprocessed files
compute_age_unc = ComputeAgeUncertainty(mesh_dir, age_dir, unc_dir, force=force)
for unc_file in tqdm(unproc_unc_files):
    status = compute_age_unc(unc_file)
