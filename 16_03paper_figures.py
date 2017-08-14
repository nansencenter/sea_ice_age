import os
import glob
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

## FIGURE 3. SIA density historgrams
factors = [2,4,8]
bins = np.arange(2,10)

sia_all = []
for factor in factors:
    ifile = '/files/sea_ice_age/fowler_nsidc_1985_f%02d/sia/1984-12-30_sia.npz' % factor
    sia = np.load(ifile)['sia']
    sia_all.append(sia[sia > 1])

plt.hist(sia_all, bins)
locs, vals = plt.yticks()
plt.yticks(locs, locs*100./1000000.)
plt.legend(['x2', 'x4', 'x8'], title='density', loc=0, ncol=3)
plt.xlabel('Sea Ice Age, years')
plt.ylabel('area, $10^6 km^2$')
plt.savefig('figure_03_sia_dens_histo.png', dpi=150, bbox_inches='tight', pad_inches=0)
plt.close()
