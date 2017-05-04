import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *
from iceagelib import *

### MAKE SIA MAPS

# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')

## FOWLER NSIDC
save_max_age('/files/sea_ice_age/fowler_nsidc_1985_f02/', 1)
save_max_age('/files/sea_ice_age/fowler_nsidc_1985_f04/', 2)
save_max_age('/files/sea_ice_age/fowler_nsidc_1985_f08/', 4)
save_max_age('/files/sea_ice_age/fowler_nsidc_2016_f02/', 1)

# source domain
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_sia_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')

prod = 'sia'

idir = '/files/sea_ice_age/fowler_nsidc_2016_f02/sia/'
ifiles = sorted(glob.glob(idir + '2015-12*.npz'))
[make_map(ifile, prod, nsidc_sia_dom, dst_dom) for ifile in ifiles]

### 
ifiles =[
'/files/sea_ice_age/fowler_nsidc_1985_f02/sia/1984-12-30_sia.npz',
'/files/sea_ice_age/fowler_nsidc_1985_f04/sia/1984-12-30_sia.npz',
'/files/sea_ice_age/fowler_nsidc_1985_f08/sia/1984-12-30_sia.npz',
'/files/sea_ice_age/fowler_nsidc_1985_f02/sia/1984-09-09_sia.npz',
'/files/sea_ice_age/fowler_nsidc_1985_f04/sia/1984-09-09_sia.npz',
'/files/sea_ice_age/fowler_nsidc_1985_f08/sia/1984-09-09_sia.npz',
]
[make_map(ifile, prod, nsidc_sia_dom, dst_dom, vmax=8) for ifile in ifiles]


## FOWLER OSI
save_max_age('/files/sea_ice_age/fowler_osi_fv1/', 1)

# source domain
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

prod = 'sia'
idir = '/files/sea_ice_age/fowler_osi_fv1/sia/'

ifiles = sorted(glob.glob(idir + '2015-12*.npz'))
[make_map(ifile, prod, osi_sic_dom, dst_dom) for ifile in ifiles]

