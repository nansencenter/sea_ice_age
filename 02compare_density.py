import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from ovl_plugins.lib.lagrangian import rungekutta4

#### IMPLEMENT THE NSIDC ICE AGE ALGORITHM
from iceagelib import *

idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
idir_ia = '/files/nsidc0611_seaice_age_v3/'
#idir = 

#icemotion.grid.week.1979.07.n.v3.bin
ifiles = sorted(glob.glob(idir_uv + '*.bin'))[:321]

# f2
ice_age_f2 = collect_age(1984, 52, ifiles, idir='/files/sea_ice_age/nsidc_f2/')

## REDUCE RESOLUTION OF NERSC PRODUCT 2 times
ice_age_f4 = collect_age(1984, 52, ifiles, idir='/files/sea_ice_age/nsidc_f4/')
ice_age_f4 = np.nanmax(np.stack([ice_age_f4[0::2, 0::2],
                                 ice_age_f4[1::2, 1::2]]), axis=0)

## REDUCE RESOLUTION OF NERSC PRODUCT 4 times
ice_age_f8 = collect_age(1984, 52, ifiles, idir='/files/sea_ice_age/nsidc_f8/')
ice_age_f8 = np.nanmax(np.stack([ice_age_f8[0::4, 0::4],
                                 ice_age_f8[1::4, 1::4],
                                 ice_age_f8[2::4, 2::4],
                                 ice_age_f8[3::4, 3::4]]), axis=0)

# get NCIDC ice age
nsidc_f = sorted(glob.glob(idir_ia + '*.bin'))[51]
nsidc_age = np.fromfile(nsidc_f, np.uint8).reshape(361*2,361*2).astype(np.float16)
nsidc_age[nsidc_age==255] = np.nan
nsidc_age /= 5.

ice_age_f2[np.isnan(nsidc_age)] = np.nan
ice_age_f4[np.isnan(nsidc_age)] = np.nan
ice_age_f8[np.isnan(nsidc_age)] = np.nan

f = plt.figure(figsize=(5,5), dpi=300)

plt.subplot(2,2,1)
plt.title('NERSC, 1_density', fontsize=5)
plt.imshow(ice_age_f2, vmin=0, vmax=8, interpolation='nearest')
plt.xticks([]);plt.yticks([])

plt.subplot(2,2,2)
plt.title('NERSC, 2_density', fontsize=5)
plt.imshow(ice_age_f4, vmin=0, vmax=8, interpolation='nearest')
plt.xticks([]);plt.yticks([])

plt.subplot(2,2,3)
plt.title('NERSC, 4_density', fontsize=5)
plt.imshow(ice_age_f8, vmin=0, vmax=8, interpolation='nearest')
plt.xticks([]);plt.yticks([])

plt.subplot(2,2,4)
plt.title('NSIDC', fontsize=5)
plt.imshow(nsidc_age, vmin=0, vmax=8, interpolation='nearest')
plt.xticks([]);plt.yticks([])
plt.tight_layout(w_pad=0)
plt.savefig('iceage_nersc_nsidc_resolutions.png', dpi=300, bbox_inches='tight', pad_inches=0)
plt.close()

h_f2 = plt.hist(ice_age_f2[np.isfinite(ice_age_f2)], range(11))[0]
h_f4 = plt.hist(ice_age_f4[np.isfinite(ice_age_f4)], range(11))[0]
h_f8 = plt.hist(ice_age_f8[np.isfinite(ice_age_f8)], range(11))[0]
h_nsidc = plt.hist(nsidc_age[np.isfinite(nsidc_age)], range(11))[0]
plt.close('all')
bins = np.arange(10)
plt.bar(bins[2:], h_f2[2:], 0.2, color='b', label='x1')
plt.bar(bins[2:]+0.2, h_f4[2:], 0.2, color='g', label='x2')
plt.bar(bins[2:]+0.4, h_f8[2:], 0.2, color='r', label='x4')
plt.bar(bins[2:]+0.6, h_nsidc[2:], 0.2, color='k', label='NSIDC')
plt.legend(loc=9)
plt.savefig('iceage_nersc_nsidc_histograms.png', dpi=75, bbox_inches='tight', pad_inches=0)
plt.close()
