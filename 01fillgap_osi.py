import glob
import os

import matplotlib.pyplot as plt
from nansat import Domain, NSR
import numpy as np

from iceagelib import fill_med_osi_uv_2, fill_gap_sic, get_osi_uvf, get_osi_date, get_osi_sic, plot_pngs

sid_dir = '/Data/sat/downloads/osi405c_demo_archive'
sic_dir = '/Data/sat/downloads/osi_ice_conc_nh_postere'

# OSISAF ICE CONC
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

# OSISAF ICE DRIFT
sid_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 62500 62500')

### MEDIAN FILTER PRIOR TO FILLING GAPS
sid_dom2x = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 31250 31250')
sid_dom4x = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 15625 15625')

w = 0.5
cMin = 15

fill_func = fill_med_osi_uv_2
odir = '/data1/antonk/sea_ice_age/osi405c_demo_archive_filled_v6/'
odir_png = f'{odir}/png'
min_conc = 0.05
n = 3
y,x = np.ogrid[-n:n+1, -n:n+1]
mask = x*x + y*y <= n*n
print(mask.shape)
print(mask)
footprint = [mask, mask]
zf = 2
sid_dom_x = sid_dom2x

# process 10 files for testing
start = 0
stop = 10
force = False
glob_mask = '201[2,3,4,5,6,7]*'

## 1. FILL GAPS IN SIC
dates2fix = [
    {'A': 20121126, 'B': 20121128, 'C': 20121127},
    {'A': 20121205, 'B': 20121207, 'C': 20121206},
    {'A': 20121224, 'B': 20121226, 'C': 20121225},
    {'A': 20130117, 'B': 20130120, 'C': 20130118},
    {'A': 20130117, 'B': 20130120, 'C': 20130119},
    {'A': 20130415, 'B': 20130417, 'C': 20130416},
    {'A': 20130812, 'B': 20130814, 'C': 20130813},        
    {'A': 20140508, 'B': 20140512, 'C': 20140509},
    {'A': 20140508, 'B': 20140512, 'C': 20140510},
    {'A': 20140508, 'B': 20140512, 'C': 20140511},
    {'A': 20140814, 'B': 20140816, 'C': 20140815},
    {'A': 20150414, 'B': 20150416, 'C': 20150415},
    {'A': 20150629, 'B': 20150701, 'C': 20150630},
    {'A': 20160404, 'B': 20160406, 'C': 20160405},
    {'A': 20160406, 'B': 20160410, 'C': 20160407},
    {'A': 20160406, 'B': 20160410, 'C': 20160408},
    {'A': 20160406, 'B': 20160410, 'C': 20160409},
    {'A': 20160731, 'B': 20160803, 'C': 20160801},
    {'A': 20160731, 'B': 20160803, 'C': 20160802},
]
for date2fix in dates2fix:
    fill_gap_sic(sic_dir, date2fix['A'], date2fix['B'], date2fix['C'])

# 2. Increase resolution of U,V in OSISAF data and save as NPZ
sid_files = sorted(glob.glob(f'{sid_dir}/ice_drift_nh_polstere-625_multi-oi_{glob_mask}.nc'))[start:stop]
uvf = None
for i, sid_file in enumerate(sid_files):
    ofile = os.path.join(odir, os.path.basename(sid_file) + '.npz')
    if os.path.exists(ofile) and not force:
        continue    

    u, v, f = get_osi_uvf(sid_file)
    sid_date = get_osi_date(sid_file)
    c = get_osi_sic(sic_dir, sid_date)
    u1, v1, c1, uvf = fill_func(u, v, c, sid_dom_x=sid_dom_x, sic_dom=sic_dom, uvf0=uvf, footprint=footprint, min_conc=min_conc, zf=zf)

    print('Save:', i, os.path.basename(ofile))
    np.savez(ofile, u=u1, v=v1, c=c1)

#3. Save U, V, C as PNG files
for sid_file in sid_files:
    ofile = os.path.join(odir, os.path.basename(sid_file) + '.npz')
    plot_pngs(ofile, odir_png)
