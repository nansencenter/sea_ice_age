### RUN ADVECTION OF ICE AGE (FROM NSIDC or OSISAF) USING OVL_PLUGINS
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *

from iceagelib import *

idir_uv = '/files/sea_ice_age/osi405c_demo_archive_filled/'

# without concentration
odir = '/files/sea_ice_age/osi_newprop_f5_zoom1/'
conc = False

#"""
# with concentration
odir = '/files/sea_ice_age/osi_newprop_f5_zoom1_conc/'
conc = True
#"""

res = 12500
factor = 5
h = 60 * 60 * 24

ifiles = sorted(glob.glob(idir_uv + 'ice_drift_nh_polstere-625_multi-oi_*.npz'))

for year in [2012, 2013, 2014, 2015]:
    i = get_i_of_file_osi(year, 10, 1, ifiles)
    propagate_from_newprop(i, ifiles, res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                get_date=get_osi_date,
                odir=odir,
                conc=conc)
    vis_ice_npz(odir + 'icemap_%d-10-01' % year)
