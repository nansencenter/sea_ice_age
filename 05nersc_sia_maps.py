import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *

from iceagelib import *
from cmocean import cm

### MAKE SIA/MYI MAPS with NERSC method (weighted average)

# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')

## NSIDC U/V
"""
nsidc_sid_dir = '/files/nsidc0116_icemotion_vectors_v3/'
reader = get_nsidc_uv
get_date = get_nsidc_date
factor = 2

# 2016
nsidc_sia_dir = '/files/sea_ice_age/nersc_nsidc_2016/'
nsidc_sid_files = sorted(glob.glob(nsidc_sid_dir + '*201[3,4,5].??.n.v3.bin'))
save_mean_age(nsidc_sid_files, nsidc_sia_dir, reader, get_date, factor=factor)

# source domain
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_sia_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')
idir = nsidc_sia_dir + 'sia/'
ifiles = sorted(glob.glob(idir + '201[3,4,5]-??-??_sia.npz'))

prod = 'sia'
for ifile in ifiles:
    make_map(ifile, prod, nsidc_sia_dom, dst_dom)

prod = 'myi'
for ifile in ifiles:
    make_map(ifile, prod, nsidc_sia_dom, dst_dom, vmin=0, vmax=1, cmap=cm.ice)
#"""

### OSI U/V
osi_sid_dir = '/files/sea_ice_age/osi405c_demo_archive_filled_v4/'
osi_sia_dir = '/files/sea_ice_age/nersc_osi_fv4_2017_conc/'
reader = get_osi_uvc_filled
get_date = get_osi_date
# source domain
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

for yy in [2016, 2017]:
    # save mean age
    osi_sid_files = sorted(glob.glob(osi_sid_dir + '*_%d*.npz' % yy))
    save_mean_age(osi_sid_files, osi_sia_dir, reader, get_date)
    # make maps (quickly)
    idir = osi_sia_dir + 'sia/'
    ifiles = sorted(glob.glob(idir + '*%d-??-??_sia.npz' % yy))
    for ifile in ifiles:
        print ifile
        plt.imsave(ifile+'_sia.png', np.load(ifile)['sia'], vmin=0, vmax=5, cmap='jet')
        plt.imsave(ifile+'_myi.png', np.load(ifile)['myi'], vmin=0, vmax=1, cmap=cm.ice)

## make maps (nicely)
#make_map(ifile, 'sia', osi_sic_dom, dst_dom, vmin=0, vmax=5, cmap='jet', text=os.path.basename(ifile).split('_')[0])
#make_map(ifile, 'myi', osi_sic_dom, dst_dom, vmin=0, vmax=1, cmap=cm.ice, text=os.path.basename(ifile).split('_')[0])
        
