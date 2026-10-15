import glob
import sys
from pathlib import Path
from tqdm import tqdm

try:
    from lmsiage.grid_age import ExportNetcdf
except ModuleNotFoundError as exc:
    if exc.name != 'lmsiage':
        raise
    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root))
    from lmsiage.grid_age import ExportNetcdf

lag_dir = 'NERSC_arctic25km_sea_ice_age_v2p1'
age_grd_dir = f'./grid'

export_config = dict(
    age_grd_dir = age_grd_dir,
    dst_root_dir = f'./cmems_nc',
    osisaf_sic_dir = '/Volumes/sim/data/OSISAF_ice_conc_CDR_v3p0',
    product_id = "SEAICE_ARC_PHY_SIAGE_L4_MY_011_025",
    file_prefix='nersc_siage_arctic65_',
    sia_name='siage',
    dataset_doi = None,
    add_met_publisher = False,
    step = 3,
    export_conc = True,
    export_uncert = False,
    export_lonlat = True,
    monthly = False,
    make_monthly_dirs = True,
    max_fraction_age = 6,
    force = False,
)

export_netcdf = ExportNetcdf(**export_config)
export_netcdf.set_attributes()

age_grd_files = sorted(glob.glob(f'{age_grd_dir}/*/grid_????????.zip'))[-100:]
print(len(age_grd_files), age_grd_files[0], age_grd_files[-1])

for age_grd_file in tqdm(age_grd_files):
    export_netcdf(age_grd_file)
