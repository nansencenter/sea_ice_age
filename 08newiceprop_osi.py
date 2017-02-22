### RUN ADVECTION OF ICE AGE (FROM NSIDC) USING OVL_PLUGINS

import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from ovl_plugins.lib.lagrangian import rungekutta4
from ovl_plugins.propagate.propagate import propagate_product_xy

from iceagelib import *

idir_uv = '/files/sea_ice_age/osi405c_demo_archive_filled/'
odir = '/files/sea_ice_age/osi_newprop_f5_zoom1/'
res = 12500
factor = 5
h = 60 * 60 * 24

ifiles = sorted(glob.glob(idir_uv + 'ice_drift_nh_polstere-625_multi-oi_*.npz'))[:900]

i = get_i_of_file_osi(2013, 9, 1, ifiles)
propagate_from_newprop(i, ifiles, res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                get_date=get_osi_date,
                odir=odir)
vis_ice_npz(odir + 'icemap_2013-09-01')

i = get_i_of_file_osi(2014, 9, 1, ifiles)
propagate_from_newprop(i, ifiles, res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                get_date=get_osi_date,
                odir=odir)
vis_ice_npz(odir + 'icemap_2014-09-01')

