### RUN ADVECTION OF ICE AGE (FROM NSIDC or OSISAF) USING OVL_PLUGINS
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *

from iceagelib import *
from scipy.ndimage.interpolation import zoom
from scipy.ndimage.filters import convolve, gaussian_filter


idir_uv = '/files/sea_ice_age/osi405c_demo_archive_filled/'

# without concentration
odir = '/files/sea_ice_age/osi_newprop_f5_zoom1/'
conc = False

#"""
# with concentration
odir = '/files/sea_ice_age/osi_newprop_f5_zoom1_conc/'
odir = '/files/sea_ice_age/osi_newprop_f10_zoom1_conc/'
conc = True
#"""

res = 6250
factor = 10
h = 60 * 60 * 24

ifiles = sorted(glob.glob(idir_uv + 'ice_drift_nh_polstere-625_multi-oi_*.npz'))
u = np.load(ifiles[1])['u']

u = fill_gaps_nn(u, 10)
u[np.isnan(u)] = 0


u0 = zoom(u, 10, order=0)
u1 = zoom(u, 10, order=1)
u2 = zoom(u, 10, order=2)

row = 920

convsize = factor
u0c = convolve(u0[row], np.ones(convsize)/convsize)
u0g = gaussian_filter(u0[row], 1)

plt.plot(u0[row], '.-', label='0')
plt.plot(u1[row], '.', label='1')
plt.plot(u2[row], '-', label='2')
plt.plot(u0c, '.-', label='0C')
plt.plot(u0g, '.-', label='0G')
plt.legend()
plt.show()


ug = gaussian_filter(u0, 1)
plt.imshow(ug);plt.show()
