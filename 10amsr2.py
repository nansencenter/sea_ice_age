import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.ndimage.filters import minimum_filter

from nansat import *

from iceagelib import *

idir = '/media/antonk/Data/files/sea_ice_age/amsr2/'
orbs = ['A', 'D']
dats = ['%02d' % d for d in range(1,11)]
ifiles = [[glob.glob('%s/GW1AM2_201506%s_01D_PNM%s_L3SGT89*.h5' % (idir, dat, orb))[0] for orb in orbs] for dat in dats]
ifiles = zip(*ifiles)[0] + zip(*ifiles)[1]

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
wm10 = minimum_filter(wm, 5)

gpi = np.isfinite(bts.sum(axis=0)) * (bts[0] > 100) * wm10

pca = np.zeros_like(bts) + np.nan
pca[:, gpi] = PCA().fit_transform(bts[:, gpi].T).T

f = Figure(pca[:3])
clim = f.clim_from_histogram()
f.process(cmin=clim[0], cmax=clim[1])
f.save('pcargb.png')

