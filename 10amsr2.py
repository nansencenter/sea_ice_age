import os
import glob

import datetime as dt

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from sklearn.decomposition import PCA
from scipy.ndimage.filters import minimum_filter

from nansat import *

from iceagelib import *


## read NSIDC SIA
idir_ia = '/files/nsidc0611_seaice_age_v3/'
nsidc_f = glob.glob(idir_ia + 'iceage.grid.week.2015.09.n.v3.bin')[0]
nsidc_age = np.fromfile(nsidc_f, np.uint8).reshape(361*2,361*2).astype(np.float32)
nsidc_age[nsidc_age==255] = np.nan
nsidc_age /= 5.
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')


idir_uv = '/files/sea_ice_age/osi405c_demo_archive_filled/'
uvfiles = sorted(glob.glob(idir_uv + 'ice_drift_nh_polstere-625_multi-oi_*.npz'))


osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'
osi_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 12500 12500')

idir = '/media/antonk/Data/files/sea_ice_age/amsr2/'
orbs = ['A']
dates = [dt.datetime(2015,3,3),
         dt.datetime(2015,3,4),
         dt.datetime(2015,3,5)]

bts_all = []
age_all = []
wm10 = None
for date in dates:
    ifile = glob.glob('%s/GW1AM2_%s_01D_PNMA_L3SGT89*.h5' % (idir, date.strftime('%Y%m%d')))[0]

    # get BT
    bts = []
    n = Nansat(ifile)
    for band in n.bands().items():
        print band[1]['name']
        bts.append(n[band[1]['name']])

    bts = np.array(bts)
    bts[bts > 655] = np.nan

   
    if wm10 is None:
        topo = Nansat('alwdgg.tif')
        topo.reproject(n, blockSize=10)
        wm = topo[1] <= 0
        wm10 = minimum_filter(wm, 10)

    # get ice age
    agew_osi, agem_osi = get_mean_age(
                         '/files/sea_ice_age/osi_newprop_f5_zoom1_conc/',
                         date)

    u,v,c = read_uv_osi_filled(uvfiles[get_i_of_file_osi(date.year, date.month, date.day, uvfiles)], factor=5)
    agew_osi, agem_osi = add_fyi(agew_osi, agem_osi, c > 0)

    c_pro = reproject_ice(osi_dom, n, c)
    agem_osi_pro = reproject_ice(osi_dom, n, agem_osi)
    agew_osi_pro = np.zeros((agew_osi.shape[0],)+agem_osi_pro.shape)
    for i in range(agew_osi.shape[0]):
        agew_osi_pro[i] = reproject_ice(osi_dom, n, agew_osi[i])

    agew_osi_pro[:,agew_osi_pro.sum(axis=0)==0] = np.nan

    bts[:, (~wm10) + (c_pro <= 0)] = np.nan
    bts_all.append(bts)

    agem_osi_pro[(~wm10) + (c_pro <= 0)] = np.nan
    age_all.append(agem_osi_pro)
    
    #nsidc_age_pro = reproject_ice(nsidc_dom, n, nsidc_age)
    #nsidc_age_pro[(~wm10) + (c_pro <= 0)] = np.nan
    #age_all.append(nsidc_age_pro)
    

bts_all = np.hstack(bts_all)
age_all = np.vstack(age_all)

bts_age = np.vstack([bts_all, age_all[None]])

gpi = np.isfinite(bts_age.sum(axis=0)) * (bts_age[14] >= 1)

#pca = np.zeros_like(bts) + np.nan
#pca[:, gpi] = PCA().fit_transform(bts[:, gpi].T).T
#f = Figure(pca[:3])
#clim = f.clim_from_histogram()
#f.process(cmin=clim[0], cmax=clim[1])
#f.save('pcargb.png')

cc = np.corrcoef(bts_age[:, gpi])
plt.imshow(cc, cmap='bwr', vmin=-1, vmax=1);plt.show()

plt.hist2d(bts_age[14, gpi], bts_age[11, gpi], 50, [[1, 5], [170, 270]], norm=LogNorm())
plt.xlabel('NERSC SIA_conc, years')
plt.ylabel('BT 89V, K')
plt.savefig('BT_89V_vs_NERSC_SIA_conc.png')
plt.close()

"""
A = np.vstack([np.ones(len(gpi[gpi])), bts[9:, gpi]]).T
B = np.linalg.lstsq(A, agem_osi_pro[gpi])[0]

A_prd = np.vstack([np.ones(1120*760), bts[9:].reshape(5, 1120*760)]).T
age_predict = np.dot(A_prd, B).reshape(agem_osi_pro.shape)
age_predict[c_pro == 0] = 0
age_predict[np.isnan(c_pro)] = np.nan
plt.imshow(age_predict);plt.show()

plt.hist2d(agem_osi_pro[gpi], age_predict[gpi], 100);plt.show()
"""
