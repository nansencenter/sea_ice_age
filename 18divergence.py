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
xOff, yOff, xSize, ySize = 40, 120, 120, 120 # twice high resolution
#xOff, yOff, xSize, ySize = 80, 240, 240, 240 # four high resolution

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
for i, sid_file in enumerate(sid_files):
    #if i != 42: continue
    date = os.path.splitext(sid_file)[0].split('-')[-1][:8]
    print date
    u, v, f = get_osi_uvf(sid_file)
    
    u = zoom_nan(u, 2) # twice high resolution
    v = zoom_nan(v, 2) # twice high resolution
    #u = zoom_nan(u, 4) # four high resolution
    #v = zoom_nan(v, 4) # four high resolution
    
    u = u[yOff:yOff+ySize, xOff:xOff+xSize]
    v = v[yOff:yOff+ySize, xOff:xOff+xSize]
    namestr = ''
    spd_all.append(np.hypot(u,v)[np.isfinite(u)])

    diver, shear, defor = get_divergence_shear(u, v)
    trn0_def = defor[:, trn]
    trn0_u = u[:, trn]
    trn0_v = v[:, trn]
        
    u = nan_median_filter(u, 7)
    v = nan_median_filter(v, 7)
    namestr = '_med'
    spd_med_all.append(np.hypot(u,v)[np.isfinite(u)])

    #u = nan_bilateral_filter(u, 3)
    #v = nan_bilateral_filter(v, 3)
    
    #u1, v1 = zoom_med(u, v)
    #u2, v2 = zoom_med(u1, v1)
    #namestr = '_med'
    
    diver, shear, defor = get_divergence_shear(u, v)

    trn1_def = defor[:, trn]
    trn1_u = u[:, trn]
    trn1_v = v[:, trn]
    x = np.arange(len(trn0_def))
    
    """
    fig, ax1 = plt.subplots()
    ax1.plot(x, trn1_u, 'k-')
    ax1.plot(x[::2], trn0_u[::2], 'r.')
    ax1.set_ylim([-0.15, 0.5])
    ax1.set_ylabel('velocity, m/s', color='r')
    ax2 = ax1.twinx()
    ax2.plot(x, trn1_def, 'k-', label='smoothed')
    ax2.plot(x[::2], trn0_def[::2], 'b.', label='original')
    ax2.set_ylim([-14, 14])
    ax2.set_ylabel('deformation, %/day', color='b')
    ax2.set_xlim([0, 80])
    plt.legend()
    #plt.show()
    #raise
    plt.savefig('transect_%0003d.png' % i)
    plt.close()
    """
    
    """
    nmap = Nansatmap(sid_n)
    nmap.imshow(defor, cmap=cm.tempo, vmin=0, vmax=10)
    nmap.quiver(u, v, step=4)
    nmap.add_colorbar(shrink=0.5, fontsize=10, format='%d')
    plt.text(2200000, 500000, date)
    nmap.save('divergence%s_%0003d.png' % (namestr, i), dpi=150)
    """

"""
for ifile in sorted(glob.glob('divergence%s_???.png' % namestr)):
    print ifile
    a = plt.imread(ifile)
    a = a[0:(a.shape[0]/2)*2, 0:(a.shape[1]/2)*2, :]
    plt.imsave(ifile, a)
#avconv -y -r 4 -i divergence_med_%0003d.png -r 25 -c:v libx264 -crf 20  -pix_fmt yuv420p -b 5000 divergence_med.mov
"""
"""
ifiles = sorted(glob.glob('/files/sea_ice_age/divergence/res62500/*.png'))
idirs = ['/files/sea_ice_age/divergence/' + wdir for wdir in
        ['res62500/', 'res62500/', 'res31250/', 'res15625_med20/']]

for ifile in ifiles:
    filename = os.path.basename(ifile)
    print filename
    call(['montage', idirs[0]+filename,
                     idirs[1]+filename.replace('_', '_med_'),
                     idirs[2]+filename.replace('_', '_med_'),
                     idirs[3]+filename.replace('_', '_med_'),
                     '-tile', '2x2', '-geometry', '+0+0', filename])
"""
spd_all = np.hstack(spd_all)
spd_med_all = np.hstack(spd_med_all)
plt.hist(np.log10(spd_all), 100, label='original')
plt.hist(np.log10(spd_med_all), 100, alpha=0.5, label='smoothed')
plt.legend()
plt.xlabel('LOG10[speed, m/s]')
plt.savefig('smoothed_speed_distribution.png', dpi=150)
plt.close('all')
