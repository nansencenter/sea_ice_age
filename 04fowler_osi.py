import glob
from nansat import *

from iceagelib import *

'''
res = 62500
h = 60 * 60 * 24
factor=5


# OLD OSI DRIFT FILLED
odir='/files/sea_ice_age/osi_f5/'

osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'
osi_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr %f %f' % (res/factor, res/factor))

nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')

ifiles = sorted(glob.glob('/files/sea_ice_age/osi405c_demo_archive_filled/*npz'))


propagate_from(get_i_of_file_osi(2012, 10, 1, ifiles), ifiles[:1000], res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                repro=(osi_dom, nsidc_dom),
                odir=odir)

propagate_from(get_i_of_file_osi(2013, 9, 15, ifiles), ifiles[:1000], res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                repro=(osi_dom, nsidc_dom),
                odir=odir)

propagate_from(get_i_of_file_osi(2014, 9, 15, ifiles), ifiles[:1000], res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                repro=(osi_dom, nsidc_dom),
                odir=odir)
'''


### NEW FILLING OF OSI SID
odir='/files/sea_ice_age/osi_newf/'

sid_dir = '/files/osi405c_demo_archive/'
sic_dir = '/files/sea_ice_age/ice_conc_grb/'

sid_files = sorted(glob.glob(sid_dir + 'ice_drift_nh_polstere-625_multi-oi_201[2,3,4,5,6,7]*.nc'))

# OSISAF Domains
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')
osi_sid_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 62500 62500')


res = 10000
h = 60 * 60 * 24
factor = 1
sigma = 2

propagate_from(get_i_of_file_osi(2012, 10, 1, sid_files), sid_files[:100], res=res,
                factor=factor, h=h, reader=get_uvc_filled_grb,
                odir=odir,
                sid_dom=osi_sid_dom, sic_dir=sic_dir, sic_dom=osi_sic_dom, sigma=sigma)
vis_ice_npz(odir + 'ice_drift_nh_polstere-625_multi-oi_201210011200', vmin=-1, vmax=1)
