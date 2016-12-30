import glob

from scipy.ndimage.filters import generic_filter
from scipy.ndimage.filters import gaussian_filter

from nansat import *

from iceagelib import *

w = 0.5
cMin = 15

ifiles = sorted(glob.glob('/Data/sat/downloads/osi405c_demo_archive/*'))

# domain of OSISAF concentration
nsr = nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
dom = Domain(nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

# read all u, v, conc
uvc = np.stack([read_uvc_osi(ifiles[i], dom) for i in range(len(ifiles))])
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
    ofile = os.path.basename(ifile).replace('.nc', '.npz')
    print 'Save', ofile
    np.savez_compressed(ofile, ice=cuv[0][i], u=cuv[1][i], v=cuv[2][i])

#vis_drift_npz('./', './')
