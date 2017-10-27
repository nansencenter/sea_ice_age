import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from ovl_plugins.lib.lagrangian import rungekutta4
from ovl_plugins.propagate.propagate import propagate_product_xy

from iceagelib import *

def get_c0_2012(osi_sic_dir):
    ofile = 'c0_2012.npy'
    if os.path.exists(ofile):
        return np.load(ofile)
    c = [get_osi_sic(osi_sic_dir, dt.datetime(2012,9,d)) for d in range(1,31)]
    c = np.array(c)
    cm = nan_median_filter(c, 5)
    c0 = np.nanmin(cm, axis=0)
    np.save(ofile, c0)
    return c0

#### RUN THE NERSC PROPAGATION
### OSI U/V
osi_sic_dir = '/files/osi_ice_conc_nh_postere/'
osi_sid_dir = '/files/sea_ice_age/osi405c_demo_archive_filled_v5/'
reader = get_osi_uvc_filled
get_date = get_osi_date

osi_sid_files = sorted(glob.glob(osi_sid_dir + '*.npz'))

## TEST
"""
odir = './'
i_start = get_osi_i_of_file(2012, 10, 1, osi_sid_files)
i_end = get_osi_i_of_file(2012, 12, 31, osi_sid_files)
propagate_nersc(i_start, i_end, osi_sid_files, reader, get_date, odir, conc=False)
vis_ice_npz(odir + 'icemap_2012')
raise
#"""

## 2012 - 2017
years = [2012, 2013, 2014, 2015, 2016]
months = [10, 9, 9, 9, 9]
days = [1, 1, 1, 1, 1]

#years = [2015, 2016]
#months = [9, 9]
#days = [15, 15]

i_end = get_osi_i_of_file(2017, 5, 31, osi_sid_files)

#conc = False
#odir = '/files/sea_ice_age/nersc_osi_fv2_2017/'

#conc = True
#odir = '/files/sea_ice_age/nersc_osi_fv4_2017_conc/'

#conc = False
#odir = '/files/sea_ice_age/nersc_osi_fv4_2017/'

#conc = True
#odir = '/files/sea_ice_age/nersc_osi_fv5_2017_conc/'

conc = False
odir = '/files/sea_ice_age/nersc_osi_fv5_2017/'

if not os.path.exists(odir):
    os.makedirs(odir)
for yy, mm, dd in zip(years, months, days):
    i_start = get_osi_i_of_file(yy, mm, dd, osi_sid_files)
    if yy == 2012:
        c0 = get_c0_2012(osi_sic_dir)
    else:
        c0 = None
    propagate_nersc(i_start, i_end, osi_sid_files, reader, get_date, odir, conc=conc, c0=c0)
    vis_ice_npz(odir + 'icemap_%04d' % yy)

