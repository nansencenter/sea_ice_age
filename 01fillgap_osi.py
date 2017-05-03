import glob
from ovl_plugins.lib.interpolation import fill_gaps_nn
from scipy import ndimage

from scipy.ndimage.filters import generic_filter
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import maximum_filter

from nansat import *

from iceagelib import *

## FILL GAPS IN OSI-SAF SEA ICE DRIFT

# get sea ice concentration
#wget -r  -nc -nd -A 'ice_conc_nh*.grb.gz' ftp://osisaf.met.no/archive/ice/conc/2017/
#wget -r  -nc -nd -w 0.5 -A 'ice_conc_nh_20*.grb.gz' ftp://osisaf.met.no/archive/ice/conc/2014/

sid_dir = '/files/osi405c_demo_archive/'
sic_dir = '/files/sea_ice_age/ice_conc_grb/'

w = 0.5
cMin = 15
odir = '/files/sea_ice_age/osi405c_demo_archive_filled_tst/'

# OSISAF ICE CONC
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

# OSISAF ICE DRIFT
sid_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 62500 62500')

start = 0
stop = None
## U/V files
sid_files = sorted(glob.glob(sid_dir + 'ice_drift_nh_polstere-625_multi-oi_201[2,3,4,5,6,7]*.nc'))[start:stop]#[1585:]

force = True
us, vs, cs = [], [], []
for i, sid_file in enumerate(sid_files):
    u, v, f = get_osi_uvf(sid_file)
    sid_date = get_osi_date(sid_file)
    c = get_osi_sic(sic_dir, sid_date)
    u, v, c = fill_osi_uv(u, v, c, sid_dom=sid_dom, sic_dom=sic_dom)
    us.append(u)
    vs.append(v)
    cs.append(c)
    
    if len(us) < 5:
        ofile = os.path.basename(sid_file) + '.npz'
        if not os.path.exists(os.path.join(odir, ofile)):
            print 'Short:', ofile
            np.savez(os.path.join(odir, ofile), u=u, v=v, c=c)
        continue
    
    ofile = os.path.basename(sid_files[i-2]) + '.npz'
    if os.path.exists(os.path.join(odir, ofile)) and not force:
        continue

    print 'Long:', ofile
    usa, vsa, csa = map(np.array, (us, vs, cs))

    # find nearest neigbour distance
    dst, ind = ndimage.distance_transform_edt(np.isnan(usa),
                                     return_distances=True,
                                     return_indices=True)

    for uv in [usa, vsa]:
        # fill with nearest neaigbour
        uvf = uv[tuple(ind)]
        # smooth
        uvf = nangaussian_filter(uvf, 3)
        # remove water
        uvf[csa == 0] = np.nan
        # remove land
        uvf[np.isnan(csa)] = np.nan
        # fill inn gaps in original U/V
        uv[np.isnan(uv)] = uvf[np.isnan(uv)]
    
    np.savez(os.path.join(odir, ofile), u=usa[2], v=vsa[2], c=csa[2])
    us.pop(0)
    vs.pop(0)
    cs.pop(0)
        
ofiles = sorted(glob.glob(odir + '*.npz'))[start:stop]
for ofile in ofiles:
    pngfile = odir + 'png/' + os.path.basename(ofile)+'_spd.png'
    if not os.path.exists(pngfile) or force:
        print 'PNG:', os.path.basename(ofile)
        u = np.load(ofile)['u']
        v = np.load(ofile)['u']
        spd = np.hypot(u, v)
        plt.imsave(pngfile, spd)
    
