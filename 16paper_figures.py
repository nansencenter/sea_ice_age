import os
import glob
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from ovl_plugins.lib.lagrangian import rungekutta4
from nansat import *
from iceagelib import *

from cmocean import cm

## FIGURE 2. COMPARISON OF SIA with different density
"""
montage\
 /files/sea_ice_age/fowler_nsidc_1985_f02/sia/1984-12-30_sia.npz_sia.png\
 /files/sea_ice_age/fowler_nsidc_1985_f04/sia/1984-12-30_sia.npz_sia.png\
 /files/sea_ice_age/fowler_nsidc_1985_f08/sia/1984-12-30_sia.npz_sia.png\
 -tile 3x1 -geometry +0+0 figure_02_sia_compar_density.png 
"""

# colorbar for SIA
cmap = cm.thermal_r
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = np.linspace(1,9,9)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

fig = plt.figure(figsize=(8, 1))
ax2 = fig.add_axes([0.05, 0.5, 0.9, 0.3])
cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i', orientation='horizontal')
cb.set_label('Sea Ice Age, years', size=12)
plt.savefig('figure_02_sia_legend.png', dpi=150, bbox_inches='tight', pad_inches=0)
plt.close()

## FIGURE 3. SIA density historgrams
factors = [2,4,8]
bins = np.arange(2,10)

sia_all = []
for factor in factors:
    ifile = '/files/sea_ice_age/fowler_nsidc_1985_f%02d/sia/1984-12-30_sia.npz' % factor
    sia = np.load(ifile)['sia']
    sia_all.append(sia[sia > 1])

plt.hist(sia_all, bins)
locs, vals = plt.yticks()
plt.yticks(locs, locs*100./1000000.)
plt.legend(['x2', 'x4', 'x8'], title='density', loc=0, ncol=3)
plt.xlabel('Sea Ice Age, years')
plt.ylabel('area, $10^6 km^2$')
plt.savefig('figure_03_sia_dens_histo.png', dpi=150, bbox_inches='tight', pad_inches=0)
plt.close()
raise

## FIGURE 4. Components of FOWLER/NSIDC


## FIGURE 8. COMPARISON OF SIA with different methods
"""
montage\
 /files/sea_ice_age/fowler_nsidc_2016_f02/sia/2015-12-31_sia.npz_sia.png\
 /files/sea_ice_age/fowler_osi_fv1/sia/2015-12-31_sia.npz_sia.png\
 /files/sea_ice_age/nersc_osi_fv1_2017/sia/2016-01-01_sia.npz_sia.png\
 /files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2016-01-01_sia.npz_sia.png\
 -tile 2x2 -geometry +0+0 figure_08_sia_compar_methods.png 
"""





"""
#### PLOT RESULTS OF THE NSIDC ICE AGE ALGORITHM

idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
idir_ia = '/files/nsidc0611_seaice_age_v3/'

# OSISAF
osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'
osi_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 6250 6250')

# NSIDC ICE AGE
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')

# DESTINATION DOMAIN    
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
#dst_dom = Domain(dst_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 6250 6250')
dst_dom = Domain(dst_nsr, '-te -2325000 -2675000 2050000 2325000 -tr 6250 6250')


# get RAW NCIDC ice age
#get_i_of_file_nsidc(1978, 44, ifiles)
nsidc_f = sorted(glob.glob(idir_ia + '*.bin'))[177]
nsidc_age = np.fromfile(nsidc_f, np.uint8).reshape(361*2,361*2).astype(np.float32)
nsidc_age[nsidc_age==255] = np.nan
nsidc_age /= 5.

plt.close('all')
nsidc_age_pro = reproject_ice(nsidc_dom, dst_dom, nsidc_age)
nmap = Nansatmap(dst_dom, resolution='l')
nmap.imshow(nsidc_age_pro, vmin=0, vmax=4, cmap='jet')
nmap.save('iceage_nsidc_raw_2015.22.png', dpi=250)
plt.close()

# GET OLD PROPAGATED NSIDC
odir = '/files/sea_ice_age/nsidc_f02_2015/'
rfiles = sorted(glob.glob(odir + '*2015.22.n.v3.bin_icemap.npz'), reverse=True)
ice_age = collect_age(rfiles)

ice_age[np.isnan(nsidc_age)] = np.nan
ice_age[nsidc_age == 0] = 0

ice_age_pro = reproject_ice(nsidc_dom, dst_dom, ice_age)

plt.close('all')
nmap = Nansatmap(dst_dom, resolution='l')
nmap.imshow(ice_age_pro, vmin=0, vmax=4, cmap='jet')
nmap.save('iceage_nsidc_nsidc_2015.22.png', dpi=250)
plt.close()

# GET OLD PROPAGATED OSI
odir = '/files/sea_ice_age/osi_f5/'
rfiles = sorted(glob.glob(odir + '*201506041200-201506061200.npz_icemap.npz'), reverse=True)
ice_age = collect_age(rfiles)

ice_age[np.isnan(nsidc_age)] = np.nan
ice_age[nsidc_age == 0] = 0

ice_age_pro = reproject_ice(nsidc_dom, dst_dom, ice_age)

plt.close('all')
nmap = Nansatmap(dst_dom, resolution='l')
nmap.imshow(ice_age_pro, vmin=0, vmax=4, cmap='jet')
nmap.save('iceage_osi_nsidc_2015.22.png', dpi=250)
plt.close()


# GET NEW PROPAGATED OSI - CONC
idir_icemap_nersc = '/files/sea_ice_age/osi_newprop_f10_zoom1_conc/'

agew_osi, agem_osi = get_mean_age(idir_icemap_nersc, dt.datetime(2015,6,2))
agem_osi_pro = reproject_ice(osi_dom, dst_dom, agem_osi)
agew_osi_pro = np.zeros((len(agew_osi),) + agem_osi_pro.shape)
for i in range(len(agew_osi)):
    agew_osi_pro[i] = reproject_ice(osi_dom, dst_dom, agew_osi[i])

ice_mask = (nsidc_age_pro > 0).astype(np.float32)
ice_mask[np.isnan(nsidc_age_pro)] = np.nan

agew_osi_pro, agem_osi_pro = add_fyi(agew_osi_pro, agem_osi_pro, ice_mask)

nmap = Nansatmap(dst_dom, resolution='l')
nmap.imshow(agem_osi_pro, vmin=0, vmax=4, cmap='jet')
nmap.save('iceage_osi_nersc_conc_2015.22.png', dpi=250)
plt.close()

# get NCIDC ice age for 2015, 22 (2015-05-01)
idir_icemap_nersc = '/files/sea_ice_age/osi_newprop_f10_zoom1/'

agew_osi, agem_osi = get_mean_age(idir_icemap_nersc, dt.datetime(2015,6,2))
agem_osi_pro = reproject_ice(osi_dom, dst_dom, agem_osi)
agew_osi_pro = np.zeros((len(agew_osi),) + agem_osi_pro.shape)
for i in range(len(agew_osi)):
    agew_osi_pro[i] = reproject_ice(osi_dom, dst_dom, agew_osi[i])

ice_mask = (nsidc_age_pro > 0).astype(np.float32)
ice_mask[np.isnan(nsidc_age_pro)] = np.nan

agew_osi_pro, agem_osi_pro = add_fyi(agew_osi_pro, agem_osi_pro, ice_mask)

nmap = Nansatmap(dst_dom, resolution='l')
nmap.imshow(agem_osi_pro, vmin=0, vmax=4, cmap='jet')
nmap.save('iceage_osi_nersc_2015.22.png', dpi=250)
plt.close()


'''
montage iceage_nsidc_raw_2015.22.png\
        iceage_osi_nsidc_2015.22.png\
        iceage_osi_nersc_2015.22.png\
        iceage_osi_nersc_conc_2015.22.png\
        -tile 2x2 -geometry +0+0 iceage_compare_mot.png 
'''


# GET NSIDC MYI/FYI maps
nsidc_dates = ['2013.15', '2014.15', '2015.15']
for nsidc_date in nsidc_dates:
    nsidc_f = sorted(glob.glob(idir_ia + '*%s*.bin' % nsidc_date))[0]
    nsidc_age = np.fromfile(nsidc_f, np.uint8).reshape(361*2,361*2).astype(np.float32)
    nsidc_age[nsidc_age==255] = np.nan
    nsidc_age /= 5.

    nsidc_age_pro = reproject_ice(nsidc_dom, dst_dom, nsidc_age)
    water_mask = (nsidc_age_pro == 0).astype(np.float)
    water_mask[water_mask != 1] = np.nan
    nmap = Nansatmap(dst_dom, resolution='l')
    nmap.imshow(nsidc_age_pro, vmin=1, vmax=2, cmap='magma')
    nmap.imshow(water_mask, vmin=10, vmax=20, cmap='bwr')
    nmap.save('iceage_nsidc_raw_%s.png' % nsidc_date, dpi=250)
    plt.close()

'''
!montage iceage_nsidc_raw_2013.15.png\
        iceage_nsidc_raw_2014.15.png\
        iceage_nsidc_raw_2015.15.png\
        -tile 3x1 -geometry +0+0 nsidc_MYI_mont_15.png 
'''


## GET MYI MAPS time for March 2012 - 2017
### OSI/SIC/NERSC
sic_dir = '/files/sea_ice_age/ice_conc_grb/'
sia_dir = '/files/sea_ice_age/osi_newprop_newf/'

osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')


idates = [
dt.datetime(2013, 3, 19),
dt.datetime(2014, 3, 19),
dt.datetime(2015, 3, 19),
dt.datetime(2016, 3, 19),
dt.datetime(2017, 3, 19),
]

for idate in idates:
    print idate
    sic = get_ice_conc_grb(sic_dir, idate)
    water_mask = (sic < 15).astype(np.float)
    water_mask[water_mask != True] = np.nan
    agew_nersc, agem_nersc = get_mean_age(sia_dir, idate)
    agew_nersc, agem_nersc = add_fyi(agew_nersc, agem_nersc, sic)
    myi_nersc = np.sum(agew_nersc[1:], axis=0)
    myi_nersc[np.isnan(sic)] = np.nan
    agem_nersc_pro = reproject_ice(osi_sic_dom, dst_dom, myi_nersc)
    water_mask_pro = reproject_ice(osi_sic_dom, dst_dom, water_mask)

    nmap = Nansatmap(dst_dom, resolution='l')
    nmap.imshow(agem_nersc_pro, vmin=0, vmax=1, cmap='magma')
    nmap.imshow(water_mask_pro, vmin=10, vmax=20, cmap='bwr')
    nmap.save('myi_nersc_%s.png' % idate.strftime('%Y%m%d'), dpi=250)
    plt.close()

'''
!montage myi_nersc_20130319.png\
         myi_nersc_20140319.png\
         myi_nersc_20150319.png\
         myi_nersc_20160319.png\
         myi_nersc_20170319.png\
         -tile 5x1 -geometry +0+0 nersc_MYI_mont.png 
'''


### GET MYI MAPS from 
"""
