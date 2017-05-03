### RUN ADVECTION OF ICE AGE (FROM NSIDC or OSISAF) USING OVL_PLUGINS
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *

from iceagelib import *


'''
## OLD FILLED UV
idir_uv = '/files/sea_ice_age/osi405c_demo_archive_filled/'

# without concentration
odir = '/files/sea_ice_age/osi_newprop_f10_zoom1/'
conc = False

# with concentration
odir = '/files/sea_ice_age/osi_newprop_f10_zoom1_conc/'
conc = True

res = 6250
factor = 10
h = 60 * 60 * 24

# new
order = 1
sigma = 0
min_flux = 0.0

ifiles = sorted(glob.glob(idir_uv + 'ice_drift_nh_polstere-625_multi-oi_*.npz'))[:1000]
for datetuple in [(2012, 10, 1), (2013, 9, 17), (2014, 9, 17), (2015, 9, 17)]:
    i = get_i_of_file_osi(datetuple[0], datetuple[1], datetuple[2], ifiles)
    propagate_from_newprop(i, ifiles, res=res,
                factor=factor, order=order, sigma=sigma, min_flux=min_flux,
                h=h, reader=read_uv_osi_filled,
                get_date=get_osi_date,
                odir=odir,
                conc=conc)
    #vis_ice_npz(odir + 'icemap_%d-10-01' % datetuple[0])
'''

### NEW FILLING OF OSI SID
odir = '/files/sea_ice_age/osi_newprop_newf/'

sid_dir = '/files/osi405c_demo_archive/'
sic_dir = '/files/sea_ice_age/ice_conc_grb/'

sid_dir = '/files/sea_ice_age/osi405c_demo_archive_filled_new/'
sid_files = sorted(glob.glob(sid_dir + 'ice_drift_nh_polstere-625_multi-oi_201[2,3,4,5,6,7]*.npz'))

# OSISAF Domains
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')
osi_sid_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 62500 62500')


res = 10000
h = 60 * 60 * 24
factor = 1
order = None
sigma = 2
min_flux = 0.0
conc = True

datetuples = [(2012, 10, 1), (2013, 9, 17), (2014, 9, 17), (2015, 9, 17), (2016, 9, 17)]
datetuples = [(2014, 9, 17), (2015, 9, 17), (2016, 9, 17)]

for datetuple in datetuples:
    i = get_i_of_file_osi(datetuple[0], datetuple[1], datetuple[2], sid_files)
    propagate_from_newprop(i, sid_files,
            res=res,
            factor=factor, order=order, sigma=sigma, min_flux=min_flux,
            h=h, reader=get_uvc_filled_grb,
            get_date=get_osi_date,
            odir=odir,
            conc=conc,
            sid_dom=osi_sid_dom, sic_dir=sic_dir, sic_dom=osi_sic_dom)
    vis_ice_npz(odir + 'icemap_%04d-%02d-%02d' % (datetuple[0], datetuple[1], datetuple[2]))
