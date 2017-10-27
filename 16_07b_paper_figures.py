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

### =====================
#####  FIGURE 8. Ice age
# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

sia_cmap = 'jet'
sic_dir = '/files/sea_ice_age/osi405c_demo_archive_filled_v5/'
sia_dir = '/files/sea_ice_age/nersc_osi_fv5_2017_conc/sia/'

c = np.load(sic_dir +'ice_drift_nh_polstere-625_multi-oi_201503011200-201503031200.nc.npz')['c']
sif = np.load(sia_dir + '2015-03-01_sia.npz')['sif']

vmax = 0.85
make_map('tmp_2015-03-01', 'sic', osi_sic_dom, dst_dom,
          vmin=0, vmax=vmax, cmap=cm.ice, title='Observed', array=c, text='A: $C_{TOT}$')

make_map('tmp_2015-03-01',
         '4yi', osi_sic_dom, dst_dom, title='from 2012', 
          vmin=0, vmax=vmax, cmap=cm.ice, array=sif[3], text='B: $C_{4YI}$')
make_map('tmp_2015-03-01',
         '3yi', osi_sic_dom, dst_dom, title='from 2013', 
          vmin=0, vmax=vmax, cmap=cm.ice, array=sif[2], text='C: $C_{3YI}$')
make_map('tmp_2015-03-01',
         '2yi', osi_sic_dom, dst_dom, title='from 2014', 
          vmin=0, vmax=vmax, cmap=cm.ice, array=sif[1], text='D: $C_{2YI}$')
make_map('tmp_2015-03-01',
         'fyi', osi_sic_dom, dst_dom, title='Remainder', 
          vmin=0, vmax=vmax, cmap=cm.ice, array=sif[0], text='E: $C_{FYI}$')
#raise
call(
'''montage\
    tmp_2015-03-01_sic.png\
    tmp_2015-03-01_4yi.png\
    tmp_2015-03-01_3yi.png\
    tmp_2015-03-01_2yi.png\
    tmp_2015-03-01_fyi.png\
    -tile 5x1 -geometry +0+0 figure_07B_averaging_sif.png'''
, shell=True)    

save_legend(cm.ice, np.linspace(0,90,19), 'Sea Ice Age Fraction, %', 'figure_07B_averaging_sif_legend.png', extend='max')

### VIZUALIZE SIA
water = reproject_ice(osi_sic_dom, dst_dom, (c <= 0.05).astype(float))

# set significant ages
## MIN, MAX, AVG AGE
sif_bin = (sif > 0.1).astype(float)
sif_ages = sif_bin * np.arange(1, len(sif_bin)+1).reshape(len(sif_bin),1,1)
sia_min = sif_ages.min(axis=0)
sia_max = sif_ages.max(axis=0)
sia_avg = np.mean(sif_ages, axis=0)
sia_rep = np.argmax(sif, axis=0).astype(float)+1

make_map('tmp_2015-03-01',
         'sia_max', osi_sic_dom, dst_dom, title='Maximum age',
         vmin=1, vmax=4, cmap=sia_cmap, array=sia_max, text='$A: MAX$', water=water)

make_map('tmp_2015-03-01',
         'sia_avg', osi_sic_dom, dst_dom, title='Average age',
         vmin=1, vmax=4, cmap=sia_cmap, array=sia_avg, text='$B: AVG$', water=water)

make_map('tmp_2015-03-01',
         'sia_rep', osi_sic_dom, dst_dom, title='Mode age',
         vmin=1, vmax=4, cmap=sia_cmap, array=sia_rep, text='$C: MOD$', water=water)

# weighted average
sia = np.load(sia_dir + '2015-03-01_sia.npz')['sia']
make_map('tmp_2015-03-01',
         'sia_wavg', osi_sic_dom, dst_dom, title='Weighted avg.', 
         vmin=1, vmax=4, cmap=sia_cmap, array=sia, text='$E: WAV$', water=water)

## INTERPOLATE SIF
def interp_sif(sif):
    sif_int = [sif[0]]
    for i in range(1, len(sif)):
        sif_int.append((sif_int[-1] + sif[i]) / 2)
        sif_int.append(sif[i])
    return np.array(sif_int)

sif_int = interp_sif(sif)
sif_int = interp_sif(sif_int)
sif_int = sif_int * sif.sum(axis=0) / sif_int.sum(axis=0)

# P50 (median)
ex_50 = np.cumsum(sif_int, axis=0)>0.5
age_layers = np.ones_like(sif_int) * np.arange(1, len(sif_int)+1).reshape(len(sif_int),1,1) /3
age_layers[(sif_int < 0.05) + np.isnan(sif_int)] = np.nan
age_layers[ex_50] = 0
sia_med = np.nanmax(age_layers, axis=0)

make_map('tmp_2015-03-01',
         'sia_med', osi_sic_dom, dst_dom, title='Median age',
         vmin=1, vmax=4, cmap=sia_cmap, array=sia_med, text='$D: P50$', water=water)


call(
"""montage\
    tmp_2015-03-01_sia_max.png\
    tmp_2015-03-01_sia_avg.png\
    tmp_2015-03-01_sia_rep.png\
    tmp_2015-03-01_sia_med.png\
    tmp_2015-03-01_sia_wavg.png\
    -tile 5x1 -geometry +0+0 figure_07B_averaging_sia.png"""
, shell=True)


save_legend(sia_cmap, np.linspace(0.75,4.25,15.), 'Sea Ice Age, years', 'figure_07B_averaging_sia_legend.png', format='%2.1f', extend='both')

"""
!montage\
    tmp_2014-03-01_sic.png tmp_2014-03-01_3yi.png tmp_2014-03-01_2yi.png tmp_2014-03-01_fyi.png /files/sea_ice_age/nersc_osi_fv4_2017_conc/sia/2014-03-01_sia.npz_sia.png\
    -tile 5x1 -geometry +0+0 figure_07B_sia_averaging.png 
"""
