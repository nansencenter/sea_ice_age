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
#####  FIGURE 7. Propagation and weighted average
# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

## FIRST COLUMN
c = np.load('/files/sea_ice_age/osi405c_demo_archive_filled_v1/ice_drift_nh_polstere-625_multi-oi_201210011200-201210031200.nc.npz')['c']
make_map('tmp_2012-10-01', 'sic', osi_sic_dom, dst_dom,
          vmin=0, vmax=100, cmap=cm.ice, title='2012-10-01', array=c, text='A: $C_{OBS}$')
make_map('tmp_2012-10-01', '2yi', osi_sic_dom, dst_dom,
          vmin=0, vmax=100, cmap=cm.ice, array=c, text='D: $C_{1+}$')
make_map('tmp_2012-10-01', '1yi', osi_sic_dom, dst_dom,
          vmin=100, vmax=200, cmap=cm.ice, array=c, text='G: $C_{0-1}$')
make_map('tmp_2012-10-01', '0yi', osi_sic_dom, dst_dom,
          vmin=100, vmax=200, cmap=cm.ice, array=c)

#c = np.load('/files/sea_ice_age/osi405c_demo_archive_filled_v1/ice_drift_nh_polstere-625_multi-oi_201303291200-201303311200.nc.npz')['c']
#make_map('tmp_2013-04-01', 'sic', osi_sic_dom, dst_dom,
#          vmin=0, vmax=100, cmap=cm.ice, title='2013-04-01', array=c, text='B: $C_{OBS}$')
#sif = np.load('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-03-29_sia.npz')['sif']
#make_map('tmp_2013-04-01',
#         '1yi', osi_sic_dom, dst_dom,
#          vmin=0, vmax=1, cmap=cm.ice, array=sif[0], text='H: $C_{FYI}$')
#make_map('tmp_2013-04-01',
#         '2yi', osi_sic_dom, dst_dom,
#          vmin=0, vmax=1, cmap=cm.ice, array=sif[1], text='E: $C_{2YI}$')

## SECOND COLUMN
c = np.load('/files/sea_ice_age/osi405c_demo_archive_filled_v1/ice_drift_nh_polstere-625_multi-oi_201309151200-201309171200.nc.npz')['c']
make_map('tmp_2013-09-15', 'sic', osi_sic_dom, dst_dom,
          vmin=0, vmax=100, cmap=cm.ice, title='2013-09-15', array=c, text='B: $C_{OBS}$')
sif = np.load('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-09-17_sia.npz')['sif']
make_map('tmp_2013-09-15',
         '3yi', osi_sic_dom, dst_dom,
          vmin=0, vmax=1, cmap=cm.ice, array=sif[2], text='E: $C_{2+}$')
make_map('tmp_2013-09-15',
         '2yi', osi_sic_dom, dst_dom,
          vmin=0, vmax=1, cmap=cm.ice, array=sif[1], text='H: $C_{1-2}$')
make_map('tmp_2013-09-15',
         '1yi', osi_sic_dom, dst_dom,
          vmin=100, vmax=200, cmap=cm.ice, array=c, text='K: $C_{0-1}$')


## THIRD COLUMN
c = np.load('/files/sea_ice_age/osi405c_demo_archive_filled_v1/ice_drift_nh_polstere-625_multi-oi_201409151200-201409171200.nc.npz')['c']
make_map('tmp_2014-09-15', 'sic', osi_sic_dom, dst_dom,
          vmin=0, vmax=100, cmap=cm.ice, title='2014-09-15', array=c, text='C: $C_{OBS}$')
sif = np.load('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2014-09-15_sia.npz')['sif']
make_map('tmp_2014-09-15',
         '4yi', osi_sic_dom, dst_dom,
          vmin=0, vmax=1, cmap=cm.ice, array=sif[3], text='F: $C_{3+}$')
make_map('tmp_2014-09-15',
         '3yi', osi_sic_dom, dst_dom,
          vmin=0, vmax=1, cmap=cm.ice, array=sif[2], text='I: $C_{2-3}$')
make_map('tmp_2014-09-15',
         '2yi', osi_sic_dom, dst_dom,
          vmin=0, vmax=1, cmap=cm.ice, array=sif[1], text='L: $C_{1-2}$')

"""
!montage\
    tmp_2012-10-01_sic.png tmp_2013-09-15_sic.png tmp_2014-09-15_sic.png\
    tmp_2012-10-01_2yi.png tmp_2013-09-15_3yi.png tmp_2014-09-15_4yi.png\
    tmp_2012-10-01_1yi.png tmp_2013-09-15_2yi.png tmp_2014-09-15_3yi.png\
    tmp_2012-10-01_0yi.png tmp_2013-09-15_1yi.png tmp_2014-09-15_2yi.png\
    -tile 3x4 -geometry +0+0 figure_07_weighted_sif.png 
"""
save_legend(cm.ice, np.linspace(0,100,20), 'Sea Ice Age Fraction, %', 'figure_07_sif_legend.png')


make_map('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2012-10-01_sia.npz',
         'sia', osi_sic_dom, dst_dom,
         vmin=1, vmax=4, cmap='jet', text='$SIA_{2012}$')

make_map('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-09-15_sia.npz',
         'sia', osi_sic_dom, dst_dom,
         vmin=1, vmax=4, cmap='jet', text='$SIA_{2013}$')

make_map('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2014-09-15_sia.npz',
         'sia', osi_sic_dom, dst_dom,
         vmin=1, vmax=4, cmap='jet', text='$SIA_{2014}$')


"""
!montage\
    /files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2012-10-01_sia.npz_sia.png\
    /files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-09-15_sia.npz_sia.png\
    /files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2014-09-15_sia.npz_sia.png\
    -tile 3x1 -geometry +0+0 figure_07_weighted_sia.png     
"""

save_legend('jet', np.linspace(0.,4.,13.), 'Sea Ice Age, years', 'figure_07_sia_legend.png', format='%2.1f')

