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
idir_uv = '/files/sea_ice_age/osi405c_demo_archive_filled/'

ice_avg = []
ice_std = []
ice_p16 = []
ice_p85 = []

mms = range(1, 13)
for mm in mms:
    print mm
    ifiles = sorted(glob.glob(idir_uv + 'ice_drift_nh_polstere-625_multi-oi_????%02d*' % mm))
    print len(ifiles)
    ice = np.array([np.load(ifile)['ice'] for ifile in ifiles])
    ice[ice < 15] = np.nan

    ice_avg.append(np.nanmean(ice))
    ice_std.append(np.nanstd(ice))
    ice_p16.append(np.nanpercentile(ice, 16))
    ice_p85.append(np.nanpercentile(ice, 86))

ice_avg = np.array(ice_avg)
ice_std = np.array(ice_std)
ice_p16 = np.array(ice_p16)
ice_p85 = np.array(ice_p85)

plt.fill_between(mms, ice_p16, ice_p85, alpha=0.5, color='gray')
plt.plot(mms, ice_avg, 'k.-')
plt.xlabel('month')
plt.ylabel('concentration')
plt.savefig('monthly_pack_sic.png')
plt.close()
