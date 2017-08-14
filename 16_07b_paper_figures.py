import os
import glob
import datetime as dt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import wget
import netCDF4 as nc4
import wget

from ovl_plugins.lib.lagrangian import rungekutta4
from nansat import *
from iceagelib import *

from cmocean import cm

### =====================
#####  FIGURE 8. Ice age
# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

sia_cmap = 'jet'

c = np.load('/files/sea_ice_age/osi405c_demo_archive_filled_v1/ice_drift_nh_polstere-625_multi-oi_201403011200-201403031200.nc.npz')['c']
make_map('tmp_2014-03-01', 'sic', osi_sic_dom, dst_dom,
          vmin=0, vmax=100, cmap=cm.ice, title='Observed', array=c, text='A: $C_{TOT}$')

sif = np.load('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2014-03-01_sia.npz')['sif']
make_map('tmp_2014-03-01',
         '3yi', osi_sic_dom, dst_dom, title='from 2012', 
          vmin=0, vmax=1, cmap=cm.ice, array=sif[2], text='B: $C_{3YI}$')
make_map('tmp_2014-03-01',
         '2yi', osi_sic_dom, dst_dom, title='from 2013', 
          vmin=0, vmax=1, cmap=cm.ice, array=sif[1], text='C: $C_{2YI}$')
make_map('tmp_2014-03-01',
         'fyi', osi_sic_dom, dst_dom, title='Remainder', 
          vmin=0, vmax=1, cmap=cm.ice, array=sif[0], text='D: $C_{FYI}$')
water = (c < 10).astype(float)
make_map('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2014-03-01_sia.npz',
         'sia', osi_sic_dom, dst_dom, title='SIA',
         vmin=1, vmax=3, cmap=sia_cmap, text='$E$', water=water)

save_legend(sia_cmap, np.linspace(1.,3.,11.), 'Sea Ice Age, years', 'figure_07B_sia_legend.png', format='%2.1f')

"""
!montage\
    tmp_2014-03-01_sic.png tmp_2014-03-01_3yi.png tmp_2014-03-01_2yi.png tmp_2014-03-01_fyi.png /files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2014-03-01_sia.npz_sia.png\
    -tile 5x1 -geometry +0+0 figure_07B_sia_averaging.png 
"""
