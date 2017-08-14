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

## FIGURE 4. Components of FOWLER/NSIDC
# source domain
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_sia_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')
# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -700000 -1400000 700000 1600000 -tr 10000 10000')

idir = '/files/sea_ice_age/fowler_nsidc_1985_f08/'
rfiles = sorted(glob.glob(idir + '*1984-12-30.npz'), reverse=True)
age_sums = []
for rfile in rfiles:
    age = np.load(rfile)['ice']
    age_sum = np.nansum(np.stack([age[0::4, 0::4],
                                  age[1::4, 1::4],
                                  age[2::4, 2::4],
                                  age[3::4, 3::4]]), axis=0)
    age_sums.append(age_sum)
age_sums = np.array(age_sums)
age_sum_tot = np.nansum(age_sums, axis=0)
for i, age_sum in enumerate(age_sums):
    age_frac = age_sum/age_sum_tot
    age_frac[~np.isfinite(age_frac)] = 0
    make_map('tmp_fraction_1984-12-30_%02d.npz' % i, 'age_frac', nsidc_sia_dom, dst_dom,
          vmin=0, vmax=0.8, cmap=cm.ice, title='%dYI' % (i+2), array=age_frac)

"""
!montage\
   tmp_fraction_1984-12-30_??.npz_age_frac.png\
   -tile 7x1 -geometry +0+0 figure_04_fowler_ice_frac.png 
"""
save_legend(cm.ice, np.linspace(0,80,20), 'Sea Ice Age Fraction, %', 'figure_04_sif_legend.png')


