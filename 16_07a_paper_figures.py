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
#####  FIGURE 7. Propagation
# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

## FIRST COLUMN
sif = np.load('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2012-10-02_sia.npz')['sif']
make_map('tmp_2012-10-01',
         '2yi', osi_sic_dom, dst_dom, title='2012-10-01', 
          vmin=0, vmax=1, cmap=cm.ice, array=sif[1], text='A: $C_{2YI}$')

sif = np.load('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-01-01_sia.npz')['sif']
make_map('tmp_2013-01-01',
         '2yi', osi_sic_dom, dst_dom, title='2013-01-01', 
          vmin=0, vmax=1, cmap=cm.ice, array=sif[1], text='B: $C_{2YI}$')

sif = np.load('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-04-04_sia.npz')['sif']
make_map('tmp_2013-04-01',
         '2yi', osi_sic_dom, dst_dom, title='2013-04-01', 
          vmin=0, vmax=1, cmap=cm.ice, array=sif[1], text='C: $C_{2YI}$')

sif = np.load('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-07-01_sia.npz')['sif']
make_map('tmp_2013-07-01',
         '2yi', osi_sic_dom, dst_dom, title='2013-07-01', 
          vmin=0, vmax=1, cmap=cm.ice, array=sif[1], text='D: $C_{2YI}$')

sif = np.load('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-10-01_sia.npz')['sif']
make_map('tmp_2013-10-01',
         '3yi', osi_sic_dom, dst_dom, title='2013-10-01', 
          vmin=0, vmax=1, cmap=cm.ice, array=sif[2], text='E: $C_{3YI}$')


"""
!montage\
    tmp_2012-10-01_2yi.png tmp_2013-01-01_2yi.png tmp_2013-04-01_2yi.png tmp_2013-07-01_2yi.png tmp_2013-10-01_3yi.png\
    -tile 5x1 -geometry +0+0 figure_07A_sif_propagation.png 
"""
#save_legend(cm.ice, np.linspace(0,100,20), 'Sea Ice Age Fraction, %', 'figure_07_sif_legend.png')
