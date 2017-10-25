import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *
from iceagelib import *

#### IMPLEMENT and RUN THE NSIDC PROPAGATION

#ftp://osisaf.met.no/archive_test/ice/drift_lr/osi405c_demo_archive/2013/01/ice_drift_nh_polstere-625_multi-oi_201301191200-201301211200.nc
#wget -w 1 -r -nc -nd -A 'ice_drift_nh_*.nc' -P /Data/sat/downloads/osi405c_demo_archive/ ftp://osisaf.met.no/archive_test/ice/drift_lr/osi405c_demo_archive/2016

### OSI U/V
osi_sid_dir = '/files/sea_ice_age/osi405c_demo_archive_filled_v4/'
reader = get_osi_uvc_filled
get_date = get_osi_date
osi_sid_files = sorted(glob.glob(osi_sid_dir + '*.npz'))


## TEST
#odir = './'
#i_start = get_osi_i_of_file(2013, 9, 15, osi_sid_files)
#i_end = get_osi_i_of_file(2013, 9, 30, osi_sid_files)
#propagate_fowler(i_start, i_end, osi_sid_files, reader, get_date, odir)
#vis_ice_npz(odir + 'icemap_2013')
#raise

odir = '/files/sea_ice_age/fowler_osi_fv4/'
if not os.path.exists(odir):
    os.makedirs(odir)

years = [2012, 2013, 2014, 2015, 2016]
months = [10, 9, 9, 9, 9]
days = [1, 15, 15, 15, 15]
i_end = get_osi_i_of_file(2017, 9, 29, osi_sid_files)

    
for yy, mm, dd in zip(years, months, days):
    i_start = get_osi_i_of_file(yy, mm, dd, osi_sid_files)
    propagate_fowler(i_start, i_end, osi_sid_files, reader, get_date, odir)
    vis_ice_npz(odir + 'icemap_%04d' % yy)
