import os
import glob
import datetime as dt
from subprocess import call

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
cmap = 'jet'
vmax = 5


## FIGURE 9. COMPARISON OF SIA with different methods
# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

# source domain
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_sia_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')

c = np.load('/files/sea_ice_age/osi405c_demo_archive_filled_v5/ice_drift_nh_polstere-625_multi-oi_201512311200-201601021200.nc.npz')['c']
water = reproject_ice(osi_sic_dom, dst_dom, (c <= 0.05).astype(float))

make_map('/files/sea_ice_age/fowler_nsidc_2016_f02/sia/2015-12-31_sia.npz',
         'sia', nsidc_sia_dom, dst_dom,
         vmin=1, vmax=vmax, cmap=cmap, text='A')

make_map('/files/sea_ice_age/fowler_osi_fv4/sia/2015-12-31_sia.npz',
         'sia', osi_sic_dom, dst_dom,
         vmin=1, vmax=vmax, cmap=cmap, text='B')

make_map('/files/sea_ice_age/nersc_osi_fv5_2017/sia/2015-12-31_sia.npz',
         'sia', osi_sic_dom, dst_dom,
         vmin=1, vmax=vmax, cmap=cmap, text='C', water=water)

make_map('/files/sea_ice_age/nersc_osi_fv6_2017_conc/sia/2015-12-31_sia.npz',
         'sia', osi_sic_dom, dst_dom,
         vmin=1, vmax=vmax, cmap=cmap, text='D', water=water)

call("""
montage\
 /files/sea_ice_age/fowler_nsidc_2016_f02/sia/2015-12-31_sia.npz_sia.png\
 /files/sea_ice_age/fowler_osi_fv4/sia/2015-12-31_sia.npz_sia.png\
 /files/sea_ice_age/nersc_osi_fv4_2017/sia/2016-01-01_sia.npz_sia.png\
 /files/sea_ice_age/nersc_osi_fv4_2017_conc/sia/2016-01-01_sia.npz_sia.png\
 -tile 2x2 -geometry +0+0 figure_08_sia_compar_methods_fv6.png 
""", shell=True)
save_legend(cmap, np.linspace(0.75,vmax+0.25,19.), 'Sea Ice Age, years', 'figure_08_sia_legend.png', format='%2.1f', extend='both')
