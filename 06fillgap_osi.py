import glob
from ovl_plugins.lib.interpolation import fill_gaps_nn
from scipy import ndimage

from scipy.ndimage.filters import generic_filter
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import maximum_filter

from nansat import *

from iceagelib import *

# get sea ice concentration
#wget -r  -nc -nd -A 'ice_conc_nh*.grb.gz' ftp://osisaf.met.no/archive/ice/conc/2017/
#wget -r  -nc -nd -w 0.5 -A 'ice_conc_nh_20*.grb.gz' ftp://osisaf.met.no/archive/ice/conc/2014/

sid_dir = '/files/osi405c_demo_archive/'
sic_dir = '/files/sea_ice_age/ice_conc_grb/'

w = 0.5
cMin = 15
odir = '/files/sea_ice_age/osi405c_demo_archive_filled_new/'

# OSISAF ICE CONC
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')
#x_pro, y_pro = sic_dom.get_geolocation_grids(dstSRS=NSR(osi_nsr))

# OSISAF ICE DRIFT
sid_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 62500 62500')
#x, y = sid_dom.get_geolocation_grids(dstSRS=NSR(osi_nsr))

start = 0
stop = None
## U/V files
sid_files = sorted(glob.glob(sid_dir + 'ice_drift_nh_polstere-625_multi-oi_201[2,3,4,5,6,7]*.nc'))[start:stop]#[1585:]

force = True
us, vs, cs = [], [], []
for i, sid_file in enumerate(sid_files):
    u, v, c = get_uvc_filled_grb(sid_file, sid_dom=sid_dom, sic_dir=sic_dir, sic_dom=sic_dom)
    us.append(u)
    vs.append(v)
    cs.append(c)
    
    if len(us) < 6:
        ofile = os.path.basename(sid_file) + '.npz'
        if not os.path.exists(os.path.join(odir, ofile)):
            print 'Short:', ofile
            np.savez(os.path.join(odir, ofile), u=u, v=v, c=c)
        continue
    
    ofile = os.path.basename(sid_files[i-2]) + '.npz'
    if os.path.exists(os.path.join(odir, ofile)) and not force:
        continue

    print 'Long:', ofile
    us.pop(0)
    vs.pop(0)
    cs.pop(0)
    
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
    
    np.savez(os.path.join(odir, ofile), u=usa[3], v=vsa[3], c=csa[3])
        

ofiles = sorted(glob.glob(odir + '*.npz'))[start:stop]
for ofile in ofiles:
    pngfile = 'new_fil_tst/' + os.path.basename(ofile)+'_spd.png'
    if not os.path.exists(pngfile) or force:
        print 'PNG:', os.path.basename(ofile)
        u = np.load(ofile)['u']
        v = np.load(ofile)['u']
        spd = np.hypot(u, v)
        plt.imsave(pngfile, spd)
    
    '''
    qstp = 1    
    plt.quiver(x[::qstp, ::qstp],
                y[::qstp, ::qstp],
                u[::qstp, ::qstp],
                v[::qstp, ::qstp],
                scale=5, color='r',
                )
    qstp = 5
    plt.quiver(x_pro[::qstp, ::qstp],
                y_pro[::qstp, ::qstp],
                uv_pro[0][::qstp, ::qstp],
                uv_pro[1][::qstp, ::qstp],
                scale=5)
    plt.xlim([-2000000, 2000000])
    plt.ylim([-2000000, 2000000])
    plt.xticks([])
    plt.yticks([])
    plt.savefig(os.path.basename(uvfile)+'.png', dpi=400, pad_inches=0, bbox_inches='tight')
    plt.close()
    '''

    

raise

c_pro = reproject_ice(osi_dom, dom, c)

# read all u, v, conc
uvf = np.stack(map(read_uvf_osi, ifiles))
u = uvc[:, 0, :, :]
v = uvc[:, 1, :, :]
c = uvc[:, 2, :, :]

cuv = [None, None, None]
k = 0
for prd in [c, u, v]:
    prdk = np.array(prd)
    if k == 0:
        prdk[prdk < 0] = np.nan
        
    prd0 = np.array(prdk)
    prdw = np.ones(prd.shape)
    
    prd0[np.isnan(prdk)] = 0
    prdw[np.isnan(prdk)] = 0
    
    prd0g = gaussian_filter(prd0, w)
    prdwg = gaussian_filter(prdw, w)
    prdg = prd0g / prdwg
    
    prdf = np.array(prdk)
    if k == 0:
        prdf[c == -101] = prdg[c == -101]
    else:
        prdf[(cuv[0] > cMin) * np.isnan(prdk)] = prdg[(cuv[0] > cMin) * np.isnan(prdk)]
        prdf[(cuv[0] > cMin) * np.isnan(prdf)] = 0

    cuv[k] = prdf
    k += 1

for i, ifile in enumerate(ifiles):
    ofile = odir + os.path.basename(ifile).replace('.nc', '.npz')
    print 'Save', ofile
    np.savez_compressed(ofile, ice=cuv[0][i], u=cuv[1][i], v=cuv[2][i])

vis_drift_npz(odir, odir)
