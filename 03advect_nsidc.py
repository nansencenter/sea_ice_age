### RUN ADVECTION OF ICE AGE (FROM NSIDC) USING OVL_PLUGINS

import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from ovl_plugins.lib.lagrangian import rungekutta4
from ovl_plugins.propagate.propagate import propagate_product_xy

from iceagelib import *

def get_ice_change(year, week, ifiles):
    i = get_i_of_file(year, week, ifiles)
    u, v, f0 = read_uv_nsidc(ifiles[i-1])
    u, v, f1 = read_uv_nsidc(ifiles[i])
    
    ice = np.zeros(f1.shape) + np.nan
    ice[(f0 > 0) * (f1 == 0)] = -1
    ice[(f0 > 0) * (f1 > 0)] = 0
    ice[(f0 == 0) * (f1 > 0)] = 1
    return ice

def propagate_product_xy_nn(x, y, x1, y1, prod0):
    ''' Propagate field of a product using new X, Y coordinates using nearest neigbour'''
    xMin, xMax = float(x.min()), float(x.max())
    yMin, yMax = float(y.min()), float(y.max())
    cols0 = np.round((x - xMin) / (xMax - xMin) * (x.shape[1] - 1)).astype(np.int16)
    rows0 = np.round((yMax - y) / (yMax - yMin) * (y.shape[0] - 1)).astype(np.int16)
    cols1 = np.round((x1 - xMin) / (xMax - xMin) * (x.shape[1] - 1)).astype(np.int16)
    rows1 = np.round((yMax - y1) / (yMax - yMin) * (y.shape[0] - 1)).astype(np.int16)
    #import ipdb; ipdb.set_trace()
    #raise

    gpi = ((cols1 >= 0) *
           (rows1 >= 0) *
           (cols1 < prod0.shape[0]) *
           (rows1 < prod0.shape[1]) *
           np.isfinite(cols1) *
           np.isfinite(rows1))
    prod1 = np.zeros(x.shape) + np.nan
    prod1[rows0[gpi], cols0[gpi]] = prod0[rows1[gpi], cols1[gpi]]

    return prod1


idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
idir_ia = '/files/nsidc0611_seaice_age_v3/'
res = 25000

#icemotion.grid.week.1979.07.n.v3.bin
ifiles = sorted(glob.glob(idir_uv + '*.bin'))

### DETECT WHEN MELTING MOSTLY STOPS AND FEEZING STARTS
#for i in range(35, 45):
#    plt.title(i)
#    plt.imshow(get_ice_change(1979, i, ifiles), vmin=-1, vmax=1)
#    plt.colorbar()
#    plt.show()
## 1979, 39

ny_week = 39

#propagate_from(1979, 39, ifiles[:100], imsave=True)
# read ice age from NSIDC
agef1 = sorted(glob.glob(idir_ia + '*.bin'))[-1]
age_nsidc = np.fromfile(agef1, np.uint8).reshape(361*2,361*2).astype(np.float16)
landmask = age_nsidc == 255

u0, v0, f0 = read_uv_nsidc(ifiles[0])
ice0 = np.zeros(f0.shape)  # water is zero age old
ice0[f0 > 0] = 2 # on 1978, week 44 all ice is in second year
ice0[landmask[::2, ::2]] = np.nan

x, y = np.meshgrid(range(u0.shape[1]), range(u0.shape[0], 0, -1))
x *= res
y *= res

x0, y0 = np.meshgrid(range(u0.shape[1]*2), range(u0.shape[0]*2, 0, -1))
x0 *= res / 2
y0 *= res / 2

h = 60 * 60 * 24 * 7 # 1 week in sec
k = 0
for i in range(0, len(ifiles)):
    year1 = int(os.path.basename(ifiles[i]).split('.')[3])
    week1 = int(os.path.basename(ifiles[i]).split('.')[4])
    print i, year1, week1
    # add age of ice
    if week1 == ny_week:
        ice0[ice0 > 0] += 1

    u1, v1, f1 = read_uv_nsidc(ifiles[i])
    u1[f1==0] = np.nan
    v1[f1==0] = np.nan    
    x1, y1 = rungekutta4(x, y, u0, v0, u1, v1, x, y, -h)
    ice1 = propagate_product_xy_nn(x, y, x1, y1, ice0)

    ice1[f1 == 0] = 0       # water
    ice1[(ice1 == 0) * (f1 > 0)] = 1 # young ice
    ice1[landmask[::2, ::2]] = np.nan # land
    #ice1[ice1 < 0.5] = 0    # wa

    
    plt.imsave('frame_nersc_nn_%03d.png' % k, ice1, vmin=0, vmax=8)
    ice0 = ice1
    u0 = u1
    v0 = v1
    k += 1

"""    
plt.imshow(ice1, vmin=0, vmax=8)
plt.colorbar()
plt.savefig('nersc_adv_19850101.png', dpi=150, pad_inches=0, bbox_inches='tight')
plt.close()
"""
