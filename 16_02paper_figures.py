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

## FIGURE 2. COMPARISON OF SIA with different density
"""
montage\
 /files/sea_ice_age/fowler_nsidc_1985_f02/sia/1984-12-30_sia.npz_sia.png\
 /files/sea_ice_age/fowler_nsidc_1985_f04/sia/1984-12-30_sia.npz_sia.png\
 /files/sea_ice_age/fowler_nsidc_1985_f08/sia/1984-12-30_sia.npz_sia.png\
 -tile 3x1 -geometry +0+0 figure_02_sia_compar_density.png 
"""
save_legend('jet', np.linspace(0,9,10), 'Sea Ice Age, years', 'figure_02_sia_legend.png', extend='both')
