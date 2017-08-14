import glob
from ovl_plugins.lib.interpolation import fill_gaps_nn
from scipy import ndimage

import cv2

from scipy.ndimage.filters import generic_filter
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.filters import median_filter
from nansat import *

from iceagelib import *


sid_dir = '/files/osi405c_demo_archive/'    
start = 700
stop = None

## U/V files
sid_files = sorted(glob.glob(sid_dir + 'ice_drift_nh_polstere-625_multi-oi_201[2,3,4,5,6,7]*.nc'))[start:stop]

u, v, f = get_osi_uvf(sid_files[13])

#u_fil = fill_gaps_nn(u, 3)
#u_bil = cv2.bilateralFilter(u_fil, 3, 3, 2)

vmin, vmax= -0.2, 0.2
"""
plt.subplot(1,3,1)
plt.imshow(u_fil, vmin=vmin, vmax=vmax, cmap='bwr')

plt.subplot(1,3,2)
plt.imshow(u_bil, vmin=vmin, vmax=vmax, cmap='bwr')
plt.subplot(1,3,3)
#"""
u_fil = fill_gaps_nn(u, 3)
plt.plot(u_fil[90,:], '.-', label='fil')
i = 0
for s0 in [1, 2, 4]:
    for s1 in [0.5, 1, 2]:
        u_fil = fill_gaps_nn(u, 3)
        u_bil = cv2.bilateralFilter(u_fil, 3, s0, s1)

        plt.plot(u_bil[90,:], '.-', label='%f %f' % (s0, s1))#, color=cm.thermal(i))
        plt.ylim([vmin, vmax])
        i += 10
plt.legend()
plt.show()
