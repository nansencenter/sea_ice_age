import glob
from ovl_plugins.lib.interpolation import fill_gaps_nn
from scipy import ndimage
from subprocess import call

import cv2

from scipy.ndimage.filters import generic_filter
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.filters import median_filter
from nansat import *

from iceagelib import *
RES = 62500.
RES = 31250.0  # twice high resolution
#RES = 15625.0  # four high resolution

def gradient2units(grd):
    return grd / RES

def defor2units(defor):
    return 100. * defor * 24. * 60. * 60.

def get_divergence_shear(u, v):
    spd = np.hypot(u, v)
    dudy, dudx = np.gradient(u)
    dvdy, dvdx = np.gradient(v)
    # original resolution 62500 m
    # time step: 1 day
    # units of deformation: % / day
    dudy, dudx, dvdy, dvdx = map(gradient2units, [dudy, dudx, dvdy, dvdx])
    
    diver = dudx + dvdy
    shear = np.hypot(dudx-dvdy, dudy+dvdx)
    defor = np.hypot(diver, shear)
    diver, shear, defor = map(defor2units, [diver, shear, defor])
    return diver, shear, defor

def nan_median_filter(arr, size):
    arr_fil = fill_gaps_nn(arr, size)
    arr_med = median_filter(arr_fil, size)
    arr_med[np.isnan(arr)] = np.nan
    return arr_med
    
def nan_bilateral_filter(arr, size=3, s0=4, s1=1):
    print 'fill_gaps_nn'
    arr_fil = fill_gaps_nn(arr, size)
    print arr_fil, 'cv2.bilateralFilter'
    arr_bil = cv2.bilateralFilter(arr_fil, size, s0, s1)
    arr_bil[np.isnan(arr)] = np.nan
    return arr_bil

def zoom_med(u, v, factor=2., size=5):
    uz = zoom_nan(u, factor)
    vz = zoom_nan(v, factor)
    uf = nan_median_filter(uz, size)
    vf = nan_median_filter(vz, size)
    uu = zoom_nan(uf, 1/factor)
    vv = zoom_nan(vf, 1/factor)
    return uu, vv

## FILL GAPS IN OSI-SAF SEA ICE DRIFT

# get sea ice concentration
#wget -r  -nc -nd -A 'ice_conc_nh*.grb.gz' ftp://osisaf.met.no/archive/ice/conc/2017/
#wget -r  -nc -nd -w 0.5 -A 'ice_conc_nh_20*.grb.gz' ftp://osisaf.met.no/archive/ice/conc/2014/

sid_dir = '/files/osi405c_demo_archive/'
sic_dir = '/files/sea_ice_age/ice_conc_grb/'

w = 0.5
cMin = 15
odir = '/files/sea_ice_age/osi405c_demo_archive_filled_v1/'

# OSISAF ICE CONC
xOff, yOff, xSize, ySize = 20, 60, 60, 60

# OSISAF ICE DRIFT
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
sid_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr %d %d' % (RES, RES))
sid_n = Nansat(domain=sid_dom, array=np.ones(sid_dom.shape()))
sid_n.crop(xOff, yOff, xSize, ySize)

start = 100
stop = 200
## U/V files
sid_files = sorted(glob.glob(sid_dir + 'ice_drift_nh_polstere-625_multi-oi_201[2,3,4,5,6,7]*.nc'))[start:stop]

trn = 60

spd_all = []
spd_med_all = []
for i in range(0, 100, 5):
    sid_file = sid_files[i]

    date = os.path.splitext(sid_file)[0].split('-')[-1][:8]
    print date
    u, v, f = get_osi_uvf(sid_file)
    u = u[yOff:yOff+ySize, xOff:xOff+xSize]
    v = v[yOff:yOff+ySize, xOff:xOff+xSize]

    uz = zoom_nan(u, 2)
    vz = zoom_nan(v, 2)

    um = nan_median_filter(uz, 7)
    vm = nan_median_filter(vz, 7)

    u2a = um[::2, ::2]
    v2a = vm[::2, ::2]

    u2b = zoom_nan(um, 0.5)
    v2b = zoom_nan(vm, 0.5)

    u2c = zoom_nan(uz, 0.5)
    v2c = zoom_nan(vz, 0.5)

    s = np.hypot(u, u)
    sa = np.hypot(u2a, u2a)
    sb = np.hypot(u2b, u2b)
    sc = np.hypot(u2c, u2c)

    vmax = 20
    plt.subplot(1,3,1)
    plt.imshow(100*(s - sa)/s, vmin=-vmax, vmax=vmax)
    plt.subplot(1,3,2)
    plt.imshow(100*(s - sb)/s, vmin=-vmax, vmax=vmax)
    plt.subplot(1,3,3)
    plt.imshow(100*(s - sc)/s, vmin=-vmax, vmax=vmax)
    plt.show()
