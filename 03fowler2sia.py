import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from ovl_plugins.lib.lagrangian import rungekutta4
from nansat import *
from iceagelib import *


#### COMPUTE SIA USING FOWLER METHOD (MAX AGE IN BIN)
idir = '/files/sea_ice_age/fowler_nsidc_1985_f02/'
ifiles = sorted(glob.glob(idir + '*.npz'))

idates = map(get_icemap_dates, ifiles)
dates0, dates1 = zip(*idates)
dst_dates = sorted(set(dates1))

dst_date = dst_dates[100]
msk = '%s*%s' % (idir, dst_date.strftime('%Y-%m-%d.npz'))
rfiles = sorted(glob.glob(msk), reverse=True)
ice_age = collect_age(rfiles)

"""
idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
idir_ia = '/files/nsidc0611_seaice_age_v3/'

# NSIDC ICE AGE
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')

# DESTINATION DOMAIN    
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
#dst_dom = Domain(dst_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 6250 6250')
dst_dom = Domain(dst_nsr, '-te -2325000 -2675000 2050000 2325000 -tr 6250 6250')


# get NCIDC ice age
nsidc_f = sorted(glob.glob(idir_ia + '*.bin'))[51]
nsidc_age = np.fromfile(nsidc_f, np.uint8).reshape(361*2,361*2).astype(np.float16)
nsidc_age[nsidc_age==255] = np.nan
nsidc_age /= 5.


bins = np.arange(2,10)
strt = 0
wdth = 0.2
ice_ages = []
factors = [2,4,8,10]
for factor in factors:
    # f2
    odir = '/files/sea_ice_age/nsidc_f%02d/' % factor
    rfiles = sorted(glob.glob(odir + '*%04d.%02d.n.v3.bin_icemap.npz' % (1984, 52)), reverse=True)
    ice_age = collect_age(rfiles)

    if factor == 4:
        ice_age = np.nanmax(np.stack([ice_age[0::2, 0::2],
                                      ice_age[1::2, 1::2]]), axis=0)
    elif factor == 8:
        ## REDUCE RESOLUTION OF NERSC PRODUCT 4 times
        ice_age = np.nanmax(np.stack([ice_age[0::4, 0::4],
                                      ice_age[1::4, 1::4],
                                      ice_age[2::4, 2::4],
                                      ice_age[3::4, 3::4]]), axis=0)
    elif factor == 10:
        ## REDUCE RESOLUTION OF NERSC PRODUCT 5 times
        ice_age = np.nanmax(np.stack([ice_age[0::5, 0::5],
                                      ice_age[1::5, 1::5],
                                      ice_age[2::5, 2::5],
                                      ice_age[3::5, 3::5],
                                      ice_age[4::5, 4::5]]), axis=0)

    ice_age[np.isnan(nsidc_age)] = np.nan
    ice_age[nsidc_age == 0] = 0

    ice_age_pro = reproject_ice(nsidc_dom, dst_dom, ice_age)

    plt.close('all')
    nmap = Nansatmap(dst_dom, resolution='l')
    nmap.imshow(ice_age_pro, vmin=0, vmax=8, cmap='jet')
    nmap.save('iceage_nsidc_f%02d.png' % factor, dpi=250)
    
    ice_ages.append(ice_age[np.isfinite(ice_age)])
    
plt.close('all')
plt.hist(ice_ages, bins)
locs, vals = plt.yticks()
plt.yticks(locs, locs*12500*12500/1000/1000/1000)
plt.legend(factors, title='factor', loc=0, ncol=2)
plt.xlabel('Sea Ice Age, years')
plt.ylabel('area, $10^6 km^2$')
plt.savefig('iceage_nsidc_histograms.png', dpi=150, bbox_inches='tight', pad_inches=0)
plt.close()

'''
montage iceage_nsidc_f02.png\
       iceage_nsidc_f04.png\
       iceage_nsidc_f08.png\
       -tile 3x1 -geometry +0+0 nsidc_age_mont00.png 
'''
"""
