import glob
from ovl_plugins.lib.interpolation import fill_gaps_nn
from scipy import ndimage
from dateutil.parser import parse


from scipy.ndimage.filters import generic_filter
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import maximum_filter

from nansat import *

from iceagelib import *

def fill_gap_sic(sic_dir, dateA, dateB, dateC):
    fileC = glob.glob(sic_dir + 'ice_conc*%s????.grb' % str(dateC))[0]
    fileC_npz = fileC + '.npz'
    if not os.path.exists(fileC_npz):
        cA = get_osi_sic(sic_dir, parse(str(dateA)), force_grb=True)
        cB = get_osi_sic(sic_dir, parse(str(dateB)), force_grb=True)
        
        cC = ((cA + cB) / 2)
        np.savez_compressed(fileC_npz, c=cC)

## FILL GAPS IN OSI-SAF SEA ICE DRIFT

# get sea ice concentration
#wget -r  -nc -nd -A 'ice_conc_nh*.grb.gz' ftp://osisaf.met.no/archive/ice/conc/2017/
#wget -r  -nc -nd -w 0.5 -A 'ice_conc_nh_20*.grb.gz' ftp://osisaf.met.no/archive/ice/conc/2014/

sid_dir = '/files/osi405c_demo_archive/'
sic_dir = '/files/osi_ice_conc_nh_postere/'

w = 0.5
cMin = 15
#fill_func = fill_med_osi_uv
#odir = '/files/sea_ice_age/osi405c_demo_archive_filled_v2/'
#med=7

#fill_func = fill_med_osi_uv
#odir = '/files/sea_ice_age/osi405c_demo_archive_filled_v3/'
#med=5

fill_func = fill_med_osi_uv_2
odir = '/files/sea_ice_age/osi405c_demo_archive_filled_v4/'
med=5

# OSISAF ICE CONC
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

# OSISAF ICE DRIFT
sid_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 62500 62500')

### MEDIAN FILTER PRIOR TO FILLING GAPS
sid_dom2x = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 31250 31250')

# process all
start = 0
stop = None
save_first = True
force_grb = False
glob_mask = '201[2,3,4,5,6,7]*'

# reprocess a gap
#start = 25
#stop = 35
#save_first = False
#force_grb = True
#glob_mask = '20160[7,8]*'

# reprocess remaining
start = 20
stop = None
save_first = False
force_grb = False
glob_mask = '20170[3,4,5,6,7,8,9]*'


# U/V files
sid_files = sorted(glob.glob(sid_dir + 'ice_drift_nh_polstere-625_multi-oi_%s.nc' % glob_mask))[start:stop]

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


#"""
uvf = None
for i, sid_file in enumerate(sid_files):
    ofile = os.path.join(odir, os.path.basename(sid_file) + '.npz')
    if os.path.exists(ofile) and not force_grb:
        continue

    u, v, f = get_osi_uvf(sid_file)
    sid_date = get_osi_date(sid_file)
    c = get_osi_sic(sic_dir, sid_date)
    u, v, c, uvf = fill_func(u, v, c, sid_dom_x=sid_dom2x, sic_dom=sic_dom, uvf0=uvf, med=med)

    if save_first:
        print 'Save:', i, os.path.basename(ofile)
        np.savez(ofile, u=u, v=v, c=c)
    save_first = True
#raise
#"""

#"""
for sid_file in sid_files:
    ofiles = sid_file.replace(sid_dir, odir)+'.npz'
    pngfile1 = odir + 'png/' + os.path.basename(ofile)+'_spd.png'
    pngfile2 = odir + 'png/sic/' + os.path.basename(ofile)+'_sic.png'
    print 'PNG:', os.path.basename(ofile)
    u = np.load(ofile)['u']
    v = np.load(ofile)['u']
    sic = np.load(ofile)['c']
    spd = np.hypot(u, v)
    plt.imsave(pngfile1, spd)
    plt.imsave(pngfile2, sic)
#"""
