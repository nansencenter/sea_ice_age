import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *
from iceagelib import *

### MAKE SIA MAPS

# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')

## FOWLER NSIDC
save_max_age('/files/sea_ice_age/fowler_nsidc_1985_f02/', 1)
save_max_age('/files/sea_ice_age/fowler_nsidc_1985_f04/', 2)
save_max_age('/files/sea_ice_age/fowler_nsidc_1985_f08/', 4)
save_max_age('/files/sea_ice_age/fowler_nsidc_2016_f02/', 1)

# source domain
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_sia_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')

prod = 'sia'
idir = '/files/sea_ice_age/fowler_nsidc_2016_f02/sia/'

ifiles = sorted(glob.glob(idir + '2015-12*.npz'))
[make_map(ifile, prod, nsidc_sia_dom, dst_dom) for ifile in ifiles]

## FOWLER OSI
save_max_age('/files/sea_ice_age/fowler_osi_fv1/', 1)

# source domain
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

prod = 'sia'
idir = '/files/sea_ice_age/fowler_osi_fv1/sia/'

ifiles = sorted(glob.glob(idir + '2015-12*.npz'))
[make_map(ifile, prod, osi_sic_dom, dst_dom) for ifile in ifiles]



"""
idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
idir_ia = '/files/nsidc0611_seaice_age_v3/'


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
