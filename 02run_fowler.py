import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *
from iceagelib import *

#### IMPLEMENT and RUN THE NSIDC PROPAGATION

#ftp://osisaf.met.no/archive_test/ice/drift_lr/osi405c_demo_archive/2013/01/ice_drift_nh_polstere-625_multi-oi_201301191200-201301211200.nc
#wget -w 1 -r -nc -nd -A 'ice_drift_nh_*.nc' -P /Data/sat/downloads/osi405c_demo_archive/ ftp://osisaf.met.no/archive_test/ice/drift_lr/osi405c_demo_archive/2016

idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
idir_ia = '/files/nsidc0611_seaice_age_v3/'

### NSIDC U/V
reader = get_nsidc_uv
get_date = get_nsidc_date
src_res = 25000 # spatial resolution of the source datasets (m)
h = 7 * 24 * 60 * 60 # temoral resolution (sec)

nsidc_sid_files = sorted(glob.glob(idir_uv + 'icemotion.grid.week*.bin'))
factors = [2, 4, 8]
odirmsk = 'nsidc_1985_f%02d/'

## 1979 - 1984
i_starts = [
    get_nsidc_i_of_file(1978, 44, nsidc_sid_files),
    get_nsidc_i_of_file(1979, 35, nsidc_sid_files),
    get_nsidc_i_of_file(1980, 35, nsidc_sid_files),
    get_nsidc_i_of_file(1981, 35, nsidc_sid_files),
    get_nsidc_i_of_file(1982, 35, nsidc_sid_files),
    get_nsidc_i_of_file(1983, 35, nsidc_sid_files),
    get_nsidc_i_of_file(1984, 35, nsidc_sid_files)
]
i_end = get_nsidc_i_of_file(1984, 52, nsidc_sid_files)

# for all zoom factors
for factor in factors:
    odir = '/files/sea_ice_age/' + odirmsk % factor
    if not os.path.exists(odir):
        os.makedirs(odir)
    # for all start years
    for i_start in i_starts:
        propagate_fowler(i_start, i_end, nsidc_sid_files, reader, get_date, src_res, h, factor, odir=odir)
