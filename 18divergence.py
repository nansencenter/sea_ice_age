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

RES = None
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
xOff, yOff, xSize, ySize = 40, 120, 120, 120 # twice high resolution
#xOff, yOff, xSize, ySize = 80, 240, 240, 240 # four high resolution

# OSISAF ICE DRIFT
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
sid_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 62500 62500')
sid_dom2x = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 31250 31250')

# OSISAF ICE CONC
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

RES = 10000

start = 100
stop = 200
step = 5
## U/V files
sid_files = sorted(glob.glob(sid_dir + 'ice_drift_nh_polstere-625_multi-oi_201[2,3,4,5,6,7]*.nc'))[start:stop:step]

trn = 300

spd_all = []
spd_med_all = []
for i, sid_file in enumerate(sid_files):
    date = os.path.splitext(sid_file)[0].split('-')[-1][:8]
    print date
    ### REGULAR WAY
    uRaw, vRaw, fRaw = get_osi_uvf(sid_file)

    u = reproject_ice(sid_dom, sic_dom, uRaw)
    v = reproject_ice(sid_dom, sic_dom, vRaw)

    u = nangaussian_filter(u, 2)
    v = nangaussian_filter(v, 2)

    diver, shear, defor = get_divergence_shear(u, v)
    trn0_def = defor[:, trn]
    trn0_u = u[:, trn]
    trn0_v = v[:, trn]

    ### APPLY MEDIAN FILTER
    # zoom2_median07_repro2x0_gaussian2
    #"""
    u = zoom_nan(uRaw, 2) # twice high resolution
    v = zoom_nan(vRaw, 2) # twice high resolution
    u = nan_median_filter(u, 7)
    v = nan_median_filter(v, 7)
    u = reproject_ice(sid_dom2x, sic_dom, u, 0)
    v = reproject_ice(sid_dom2x, sic_dom, v, 0)
    u = nangaussian_filter(u, 2)
    v = nangaussian_filter(v, 2)
    #"""

    # zoom2_median07_repro2x1
    """
    u = zoom_nan(uRaw, 2)
    v = zoom_nan(vRaw, 2)
    u = nan_median_filter(u, 7)
    v = nan_median_filter(v, 7)
    u = reproject_ice(sid_dom2x, sic_dom, u, 1)
    v = reproject_ice(sid_dom2x, sic_dom, v, 1)
    #"""

    diver, shear, defor = get_divergence_shear(u, v)
    trn1_def = defor[:, trn]
    trn1_u = u[:, trn]
    trn1_v = v[:, trn]
    x = np.arange(len(trn0_def))
    
    """
    fig, ax1 = plt.subplots()
    ax1.plot(x, trn1_u, 'k-')
    ax1.plot(x[::2], trn0_u[::2], 'r.')
    ax1.set_ylim([-0.2, 0.5])
    ax1.set_ylabel('velocity, m/s', color='r')
    ax2 = ax1.twinx()
    ax2.plot(x, trn1_def, 'k-', label='smoothed')
    ax2.plot(x[::2], trn0_def[::2], 'b.', label='original')
    ax2.set_ylim([-14, 14])
    ax2.set_ylabel('deformation, %/day', color='b')
    ax2.set_xlim([400, 600])
    plt.legend()
    #plt.show()
    plt.savefig('transect_%0003d.png' % i)
    plt.close()
    #"""
    
    #"""
    nmap = Nansatmap(sic_dom)
    nmap.imshow(defor, cmap=cm.tempo, vmin=0, vmax=10)
    nmap.quiver(u, v, step=10)
    #nmap.add_colorbar(shrink=0.5, fontsize=10, format='%d')
    #plt.text(2200000, 500000, date)
    nmap.save('divergence_%0003d.png' % i, dpi=300)
    #"""
