import glob
from nansat import *

from iceagelib import *

res = 62500
h = 60 * 60 * 24
factor=5

osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'
osi_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr %f %f' % (res/factor, res/factor))

nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')

ifiles = sorted(glob.glob('/files/sea_ice_age/osi405c_demo_archive_filled/*npz'))


propagate_from(get_i_of_file_osi(2013, 9, 15, ifiles), ifiles[:255], res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                repro=(osi_dom, nsidc_dom),
                odir='/files/sea_ice_age/20130915_npz/')

raise
"""
propagate_from(get_i_of_file_osi(2013, 9, 15, ifiles), ifiles[:255], res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                repro=(osi_dom, nsidc_dom),
                odir='/files/sea_ice_age/20130915_npz/')

propagate_from(get_i_of_file_osi(2014, 9, 15, ifiles), ifiles, res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                repro=(osi_dom, nsidc_dom),
                odir='/files/sea_ice_age/20140915_npz/')

propagate_from(get_i_of_file_osi(2015, 9, 15, ifiles), ifiles, res=res,
                factor=factor, h=h, reader=read_uv_osi_filled,
                repro=(osi_dom, nsidc_dom),
                odir='/files/sea_ice_age/20150915_npz/')

vis_ice_npz('/files/sea_ice_age/20130915_npz/')
vis_ice_npz('/files/sea_ice_age/20140915_npz/')
vis_ice_npz('/files/sea_ice_age/20150915_npz/')
"""

rdate = '20150115'
rfiles = (glob.glob('/files/sea_ice_age/20130915_npz/*%s1200.npz_icemap.npy.npz'%rdate) +
          glob.glob('/files/sea_ice_age/20140915_npz/*%s1200.npz_icemap.npy.npz'%rdate))

rfile = glob.glob('/files/sea_ice_age/osi405c_demo_archive_filled/*%s1200.npz'%rdate)[0]

ice_age = collect_age(rfiles)

# get C for dst date
u,v,c = read_uv_osi_filled(rfile)
factor=1
osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'
osi_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr %f %f' % (res/factor, res/factor))
cpro = reproject_ice(osi_dom, nsidc_dom, c)

# mask water
ice_age[cpro < 15] = 0
ice_age[np.isnan(cpro)] = np.nan

plt.imshow(ice_age[100:600, 100:600])
cbar = plt.colorbar()
plt.savefig('ice_age_20150115_fowler_osisaf.png', dpi=150)
plt.close()

#factor=2
#idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
#ifiles = sorted(glob.glob(idir_uv + '*201[3,4,5]*.bin'))
#propagate_from(get_i_of_file_nsidc(2013, 38, ifiles), ifiles, factor=factor)
#propagate_from(get_i_of_file_nsidc(2014, 38, ifiles), ifiles, factor=factor)
