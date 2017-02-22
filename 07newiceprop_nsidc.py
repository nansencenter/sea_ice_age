### RUN ADVECTION OF ICE AGE (FROM NSIDC) USING OVL_PLUGINS

import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from ovl_plugins.lib.lagrangian import rungekutta4
from ovl_plugins.propagate.propagate import propagate_product_xy

from iceagelib import *

idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
idir_ia = '/files/nsidc0611_seaice_age_v3/'

res = 12500
factor = 2

"""
# 1978 - 1984
# for comparison with NSIDC ice age product
odir = '/files/sea_ice_age/nsidc_f2_newprop_1985/'
ifiles = sorted(glob.glob(idir_uv + 'icemotion.grid.week.19*.bin'))

propagate_from_newprop(get_i_of_file_nsidc(1978, 44, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop_tmp/icemap_1978-11-05')

propagate_from_newprop(get_i_of_file_nsidc(1979, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop_tmp/icemap_1979-09-03')

propagate_from_newprop(get_i_of_file_nsidc(1980, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop_tmp/icemap_1980-09-02')

propagate_from_newprop(get_i_of_file_nsidc(1981, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop_tmp/icemap_1981-09-03')

propagate_from_newprop(get_i_of_file_nsidc(1982, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop_tmp/icemap_1982-09-03')

propagate_from_newprop(get_i_of_file_nsidc(1983, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop_tmp/icemap_1983-09-03')

propagate_from_newprop(get_i_of_file_nsidc(1984, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop_tmp/icemap_1979-09-02')
"""

# 2013 - 2016
# for comparison with OSI
odir = '/files/sea_ice_age/nsidc_f2_newprop_2015/'
ifiles = sorted(glob.glob(idir_uv + 'icemotion.grid.week.20*.bin'))

propagate_from_newprop(get_i_of_file_nsidc(2013, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz(odir + 'icemap_2013-09-03')

propagate_from_newprop(get_i_of_file_nsidc(2014, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz(odir + 'icemap_2014-09-03')

