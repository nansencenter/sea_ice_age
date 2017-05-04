import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from ovl_plugins.lib.lagrangian import rungekutta4
from ovl_plugins.propagate.propagate import propagate_product_xy

from iceagelib import *

#### RUN THE NERSC PROPAGATION
### NSIDC U/V
nsidc_sid_dir = '/files/nsidc0116_icemotion_vectors_v3/'

reader = get_nsidc_uv
get_date = get_nsidc_date
factor = 2

nsidc_sid_files = sorted(glob.glob(nsidc_sid_dir + 'icemotion.grid.week*.bin'))

## TEST
#odir = '/files/sea_ice_age/nersc_nsidc_test/'
#os.makedirs(odir)
#i_start = get_nsidc_i_of_file(1978, 44, nsidc_sid_files)
#i_end = get_nsidc_i_of_file(1979, 52, nsidc_sid_files)
#propagate_nersc(i_start, i_end, nsidc_sid_files, reader, get_date, odir, factor=factor)
#vis_ice_npz(odir + 'icemap_1978')
#raise

## 2012 - 2015
odir = '/files/sea_ice_age/nersc_nsidc_2016/'
years = [2012, 2013, 2014, 2015]
weeks = [35,   35,   35,   35]
i_end = get_nsidc_i_of_file(2015, 52, nsidc_sid_files)

for yy, ww in zip(years, weeks):
    if not os.path.exists(odir):
        os.makedirs(odir)
    i_start = get_nsidc_i_of_file(yy, ww, nsidc_sid_files)
    propagate_nersc(i_start, i_end, nsidc_sid_files, reader, get_date, odir, factor=factor)
    vis_ice_npz(odir + 'icemap_%04d' % yy)


