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
osi_sid_dir = '/files/sea_ice_age/osi405c_demo_archive_filled_v1/'
reader = get_osi_uvc_filled
get_date = get_osi_date
src_res = 10000 # spatial resolution of the source datasets (m)
h = 24 * 60 * 60 # temoral resolution (sec)
factor = 1

odir = '/files/sea_ice_age/fowler_osi_fv1/'
if not os.path.exists(odir):
    os.makedirs(odir)


osi_sid_files = sorted(glob.glob(osi_sid_dir + '*.npz'))

u, v, c = get_osi_uvc_filled(osi_sid_files[0], src_res, h, factor)

i_start = get_osi_i_of_file(2012, 10, 1, osi_sid_files)
i_end = get_osi_i_of_file(2012, 12, 31, osi_sid_files)
propagate_fowler(i_start, i_end, osi_sid_files, reader, get_date, src_res, h, factor, odir=odir)
vis_ice_npz(odir + 'icemap_2012')
