import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from ovl_plugins.lib.lagrangian import rungekutta4

#### IMPLEMENT THE NSIDC ICE AGE ALGORITHM
from iceagelib import *

idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
idir_ia = '/files/nsidc0611_seaice_age_v3/'


#icemotion.grid.week.1979.07.n.v3.bin
ifiles = sorted(glob.glob(idir_uv + '*.bin'))[:321]

#"""
factor=10
propagate_from(get_i_of_file_nsidc(1978, 44, ifiles), ifiles, factor=factor)
propagate_from(get_i_of_file_nsidc(1979, 35, ifiles), ifiles, factor=factor)
propagate_from(get_i_of_file_nsidc(1980, 35, ifiles), ifiles, factor=factor)
propagate_from(get_i_of_file_nsidc(1981, 35, ifiles), ifiles, factor=factor)
propagate_from(get_i_of_file_nsidc(1982, 35, ifiles), ifiles, factor=factor)
propagate_from(get_i_of_file_nsidc(1983, 35, ifiles), ifiles, factor=factor)
propagate_from(get_i_of_file_nsidc(1984, 35, ifiles), ifiles, factor=factor)
#"""

idir = '/files/sea_ice_age/nsidc_f10/'
ifiles = sorted(glob.glob(idir + '*%04d.%02d.n.v3.bin_icemap.npy.npz' % (1984, 52)), reverse=True)
ice_age = collect_age(ifiles)

## REDUCE RESOLUTION OF NERSC PRODUCT TWICE
#ice_age = np.nanmax(np.stack([ice_age[::2, ::2],
#                              ice_age[1::2, 1::2]]), axis=0)

## REDUCE RESOLUTION OF NERSC PRODUCT 4 times
#ice_age = np.nanmax(np.stack([ice_age[0::4, 0::4],
#                              ice_age[1::4, 1::4],
#                              ice_age[2::4, 2::4],
#                              ice_age[3::4, 3::4]]), axis=0)

## REDUCE RESOLUTION OF NERSC PRODUCT 5 times
ice_age = np.nanmax(np.stack([ice_age[0::5, 0::5],
                              ice_age[1::5, 1::5],
                              ice_age[2::5, 2::5],
                              ice_age[3::5, 3::5],
                              ice_age[4::5, 4::5]]), axis=0)

# get NCIDC ice age
nsidc_f = sorted(glob.glob(idir_ia + '*.bin'))[51]
nsidc_age = np.fromfile(nsidc_f, np.uint8).reshape(361*2,361*2).astype(np.float16)
nsidc_age[nsidc_age==255] = np.nan
nsidc_age /= 5.

ice_age[np.isnan(nsidc_age)] = np.nan

plt.subplot(1,2,1)
plt.title('Ice age, 1 Jan 1985, NERSC')
plt.imshow(ice_age, vmin=0, vmax=8, interpolation='nearest')
plt.xticks([]);plt.yticks([])
plt.subplot(1,2,2)
plt.title('Ice age, 1 Jan 1985, NSIDC')
plt.imshow(nsidc_age, vmin=0, vmax=8, interpolation='nearest')
#plt.colorbar(shrink=0.5)
plt.xticks([]);plt.yticks([])
plt.savefig('iceage_nersc_nsidc.png', dpi=300, bbox_inches='tight', pad_inches=0)
plt.close()
