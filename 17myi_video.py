import os
import glob
import shutil
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


### Video Comparison of MYI from NERSC and BREMEN


####============== BREMEN
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
brem_dom = Domain(osi_nsr, '-te -3843750 -5343750 3756250 5856250 -tr 12500 12500')

idir = '/files/sea_ice_age/bremensit/'
osi_sia_dir = '/files/sea_ice_age/nersc_osi_fv1_2017_conc/'

#ifiles = sorted(glob.glob(idir + 'MYI-NMYI-CORRECTION-*.nc'))
#for ifile in ifiles:
#    print ifile
#    brem_myi = nc4.Dataset(ifile).variables['NMYI'][:].reshape(896, 608)
#    date = os.path.splitext(ifile)[0].split('-')[3]
#    make_map(idir + 'tmp_%s' % date, 'bremen_myi', brem_dom, dst_dom,
#             array=brem_myi, vmin=0, vmax=100, cmap=cm.ice)

ifiles = sorted(glob.glob(idir + 'tmp_*_bremen_myi.png'))
for ifile in ifiles:
    ofile = ifile + '_join.png'
    if os.path.exists(ofile):
        continue
    print ifile
    
    datestr = os.path.basename(ifile).split('_')[1]

    sia_map_files = glob.glob(osi_sia_dir + 'sia/%s-%s-%s_sia.npz_myi.png' % (datestr[0:4], datestr[4:6], datestr[6:8]))
    if len(sia_map_files) >= 1:
        brem_map = plt.imread(ifile)
        sia_map = plt.imread(sia_map_files[0])
        join_map = np.zeros((max(brem_map.shape[0], sia_map.shape[0]) + 1,
                             brem_map.shape[1] + sia_map.shape[1] + 4,
                             4)) + np.nan
        join_map[0:sia_map.shape[0], 0:sia_map.shape[1], :] = sia_map
        join_map[0:brem_map.shape[0],
                (sia_map.shape[1]+4):(sia_map.shape[1]+4+brem_map.shape[1]),
                :] = brem_map

    plt.imsave(ofile, join_map)


ifiles = sorted(glob.glob(idir + '*_join.png'))
for i, ifile in enumerate(ifiles):
    shutil.copy(ifile, idir + 'joined_frame_%002d.png' % i)
    
#avconv -y -r 24 -i /files/sea_ice_age/bremensit/joined_frame_%002d.png -r 25 -c:v libx264 -crf 20  -pix_fmt yuv420p -b 5000 myi_compar.mov


#### MOVIE OF 2012 propagation only
ifiles = sorted(glob.glob('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/201[2,3]*_sia.npz_myi.png'))
for i, ifile in enumerate(ifiles):
    print i, ifile
    a = plt.imread(ifile)
    a = a[0:(a.shape[0]/2)*2, 0:(a.shape[1]/2)*2, :]
    plt.imsave('myi_frame_%03d.png' % i, a)

    
#avconv -y -r 24 -i myi_frame_%02d.png -r 25 -c:v libx264 -crf 20  -pix_fmt yuv420p -b 5000 myi_video.mov
