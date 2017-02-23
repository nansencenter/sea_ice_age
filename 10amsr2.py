import os
import glob

import datetime as dt

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.ndimage.filters import minimum_filter

from nansat import *

from iceagelib import *

idir_uv = '/files/sea_ice_age/osi405c_demo_archive_filled/'
uvfiles = sorted(glob.glob(idir_uv + 'ice_drift_nh_polstere-625_multi-oi_*.npz'))


osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'
osi_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 12500 12500')

idir = '/media/antonk/Data/files/sea_ice_age/amsr2/'
orbs = ['A']
dats = ['%02d' % d for d in range(1,11)]
ifiles = [glob.glob('%s/GW1AM2_20150101_01D_PNM%s_L3SGT89*.h5' % (idir, orb))[0] for orb in orbs]


# get BT
bts = []
for ifile in ifiles:
    print ifile
    n = Nansat(ifile)
    for band in n.bands().items():
        print band[1]['name']
        bts.append(n[band[1]['name']])

bts = np.array(bts)
bts[bts > 655] = np.nan

topo = Nansat('alwdgg.tif')
topo.reproject(n, blockSize=10)
wm = topo[1] <= 0
wm10 = minimum_filter(wm, 10)

# get ice age
thedate = dt.datetime(2015,1,1)
agew_osi, agem_osi = get_mean_age(
                     '/files/sea_ice_age/osi_newprop_f5_zoom1_conc/',
                     thedate)

u,v,c = read_uv_osi_filled(uvfiles[get_i_of_file_osi(2015, 1, 1, uvfiles)], factor=5)
agew_osi, agem_osi = add_fyi(agew_osi, agem_osi, c > 0)

agem_osi_pro = reproject_ice(osi_dom, n, agem_osi)
agew_osi_pro = np.zeros((agew_osi.shape[0],)+agem_osi_pro.shape)
for i in range(agew_osi.shape[0]):
    agew_osi_pro[i] = reproject_ice(osi_dom, n, agew_osi[i])

agew_osi_pro[:,agew_osi_pro.sum(axis=0)==0] = np.nan

bts_age = np.vstack([bts, agem_osi_pro.reshape((1,)+agem_osi_pro.shape), agew_osi_pro])

gpi = np.isfinite(bts_age.sum(axis=0)) * (bts_age[0] > 100) * wm10# * (agem_osi_pro > 1.1)

pca = np.zeros_like(bts) + np.nan
pca[:, gpi] = PCA().fit_transform(bts[:, gpi].T).T

f = Figure(pca[:3])
clim = f.clim_from_histogram()
f.process(cmin=clim[0], cmax=clim[1])
f.save('pcargb.png')

cc = np.corrcoef(bts_age[:, gpi])
plt.imshow(cc, cmap='bwr', vmin=-1, vmax=1);plt.show()
