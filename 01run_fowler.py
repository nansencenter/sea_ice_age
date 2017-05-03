import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from ovl_plugins.lib.lagrangian import rungekutta4

#### IMPLEMENT and RUN THE NSIDC PROPAGATION

#ftp://osisaf.met.no/archive_test/ice/drift_lr/osi405c_demo_archive/2013/01/ice_drift_nh_polstere-625_multi-oi_201301191200-201301211200.nc
#wget -w 1 -r -nc -nd -A 'ice_drift_nh_*.nc' -P /Data/sat/downloads/osi405c_demo_archive/ ftp://osisaf.met.no/archive_test/ice/drift_lr/osi405c_demo_archive/2016


from iceagelib import *

idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
idir_ia = '/files/nsidc0611_seaice_age_v3/'

#FOR 1984
#ifiles = sorted(glob.glob(idir_uv + 'icemotion.grid.week*.bin'))[:321]
#factors = [2, 4, 8,10]
#odirmsk = 'nsidc_f%02d/'

#FOR 2015
ifiles = sorted(glob.glob(idir_uv + 'icemotion.grid.week*.bin'))
factors = [2]
odirmsk = 'nsidc_f%02d_2015/'

for factor in factors:
    odir = '/files/sea_ice_age/' + odirmsk % factor

    '''
    # FOR 1984
    propagate_from(get_i_of_file_nsidc(1978, 44, ifiles), ifiles, factor=factor, odir=odir)
    propagate_from(get_i_of_file_nsidc(1979, 35, ifiles), ifiles, factor=factor, odir=odir)
    propagate_from(get_i_of_file_nsidc(1980, 35, ifiles), ifiles, factor=factor, odir=odir)
    propagate_from(get_i_of_file_nsidc(1981, 35, ifiles), ifiles, factor=factor, odir=odir)
    propagate_from(get_i_of_file_nsidc(1982, 35, ifiles), ifiles, factor=factor, odir=odir)
    propagate_from(get_i_of_file_nsidc(1983, 35, ifiles), ifiles, factor=factor, odir=odir)
    propagate_from(get_i_of_file_nsidc(1984, 35, ifiles), ifiles, factor=factor, odir=odir)
    #'''

    #'''
    # FOR 2015
    propagate_from(get_i_of_file_nsidc(2012, 37, ifiles), ifiles, factor=factor, odir=odir)
    propagate_from(get_i_of_file_nsidc(2013, 37, ifiles), ifiles, factor=factor, odir=odir)
    propagate_from(get_i_of_file_nsidc(2014, 37, ifiles), ifiles, factor=factor, odir=odir)
    #'''
