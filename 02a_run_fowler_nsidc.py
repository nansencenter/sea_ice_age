import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *
from iceagelib import *

#### RUN THE NSIDC PROPAGATION

#ftp://osisaf.met.no/archive_test/ice/drift_lr/osi405c_demo_archive/2013/01/ice_drift_nh_polstere-625_multi-oi_201301191200-201301211200.nc
#wget -w 1 -r -nc -nd -A 'ice_drift_nh_*.nc' -P /Data/sat/downloads/osi405c_demo_archive/ ftp://osisaf.met.no/archive_test/ice/drift_lr/osi405c_demo_archive/2016
#wget -nc -nd --http-user=korosov --http-password= --load-cookies mycookies.txt --save-cookies mycookies.txt --keep-session-cookies --no-check-certificate --auth-no-challenge -r -np -e robots=off https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0611_seaice_age_v3/data/2015/

### NSIDC U/V
nsidc_sid_dir = '/files/nsidc0116_icemotion_vectors_v3/'
reader = get_nsidc_uv
get_date = get_nsidc_date
nsidc_sid_files = sorted(glob.glob(nsidc_sid_dir + 'icemotion.grid.week*.bin'))

## TEST
#factor = 2
#odir = './'
#i_start = get_nsidc_i_of_file(1978, 44, nsidc_sid_files)
#i_end = get_nsidc_i_of_file(1979, 10, nsidc_sid_files)
#propagate_fowler(i_start, i_end, nsidc_sid_files, reader, get_date, odir, factor=factor)
#vis_ice_npz(odir + 'icemap_1978')
#raise

## 1979 - 1984
factors = [2, 4, 8]
odirmsk = 'fowler_nsidc_1985_f%02d/'
years = [1978, 1979, 1980, 1981, 1982, 1983, 1984]
weeks = [44, 35, 35, 35, 35, 35, 35]
i_end = get_nsidc_i_of_file(1984, 52, nsidc_sid_files)

## 2012 - 2015
factors = [2]
odirmsk = 'fowler_nsidc_2016_f%02d/'
years = [2012, 2013, 2014, 2015]
weeks = [35,   35,   35,   35]
i_end = get_nsidc_i_of_file(2015, 52, nsidc_sid_files)

# for all zoom factors
for factor in factors:
    odir = '/files/sea_ice_age/' + odirmsk % factor
    if not os.path.exists(odir):
        os.makedirs(odir)
    # for all start years
    for yy, ww in zip(years, weeks):
        i_start = get_nsidc_i_of_file(yy, ww, nsidc_sid_files)
        propagate_fowler(i_start, i_end, nsidc_sid_files, reader, get_date, odir, factor=factor)
        vis_ice_npz(odir + 'icemap_%04d' % yy)
