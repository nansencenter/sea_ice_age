import glob
from nansat import *

from iceagelib import *


res = 62500
h = 60 * 60 * 24
factor=5

ifiles = sorted(glob.glob('/Data/sat/downloads/osi405c_demo_archive/*_201[3,4,5]*'))

n0 = Nansat(ifiles[0])
osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'
osi_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr %f %f' % (res/factor, res/factor))

nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')


ifiles = sorted(glob.glob('/files/sea_ice_age/osi405c_demo_archive_filled/*npz'))
"""
propagate_from(get_i_of_file_osi(2013, 9, 15, ifiles), ifiles, res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                repro=(osi_dom, nsidc_dom),
                odir='./201309151200-201309171200_npz_fill')
"""

propagate_from(get_i_of_file_osi(2014, 9, 15, ifiles), ifiles, res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                repro=(osi_dom, nsidc_dom),
                odir='./201409151200-201409171200_npz_fill')

propagate_from(get_i_of_file_osi(2015, 9, 15, ifiles), ifiles, res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                repro=(osi_dom, nsidc_dom),
                odir='./201509151200-201509171200_npz_fill')

vis_ice_npz('201309151200-201309171200_npz_fill/')


factor=2
idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
ifiles = sorted(glob.glob(idir_uv + '*201[3,4,5]*.bin'))
propagate_from(get_i_of_file_nsidc(2013, 38, ifiles), ifiles, factor=factor)
propagate_from(get_i_of_file_nsidc(2014, 38, ifiles), ifiles, factor=factor)
