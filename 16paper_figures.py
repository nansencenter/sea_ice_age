import os
import glob
import datetime as dt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import wget
import netCDF4 as nc4
import wget

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
save_legend(cm.thermal_r, np.linspace(1,9,9), 'Sea Ice Age, years', 'figure_02_sia_legend.png')


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
# source domain
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_sia_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')
# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -700000 -1400000 700000 1600000 -tr 10000 10000')

idir = '/files/sea_ice_age/fowler_nsidc_1985_f08/'
rfiles = sorted(glob.glob(idir + '*1984-12-30.npz'), reverse=True)
age_sums = []
for rfile in rfiles:
    age = np.load(rfile)['ice']
    age_sum = np.nansum(np.stack([age[0::4, 0::4],
                                  age[1::4, 1::4],
                                  age[2::4, 2::4],
                                  age[3::4, 3::4]]), axis=0)
    age_sums.append(age_sum)
age_sums = np.array(age_sums)
age_sum_tot = np.nansum(age_sums, axis=0)
for i, age_sum in enumerate(age_sums):
    age_frac = age_sum/age_sum_tot
    age_frac[~np.isfinite(age_frac)] = 0
    make_map('tmp_fraction_1984-12-30_%02d.npz' % i, 'age_frac', nsidc_sia_dom, dst_dom,
          vmin=0, vmax=0.8, cmap=cm.ice, title='%dYI' % (i+2), array=age_frac)

"""
!montage\
   tmp_fraction_1984-12-30_??.npz_age_frac.png\
   -tile 7x1 -geometry +0+0 figure_04_fowler_ice_frac.png 
"""
save_legend(cm.ice, np.linspace(0,80,20), 'Sea Ice Age Fraction, %', 'figure_04_sif_legend.png')


#### FIGURE 5. AVERAGE SIC 
idir_uv = '/files/sea_ice_age/osi405c_demo_archive_filled_v1/'
ice_avg = []
ice_std = []
ice_p16 = []
ice_p85 = []

mms = range(1, 13)
for mm in mms:
    print mm
    ifiles = sorted(glob.glob(idir_uv + 'ice_drift_nh_polstere-625_multi-oi_????%02d*' % mm))
    print len(ifiles)
    ice = np.array([np.load(ifile)['c'] for ifile in ifiles])
    ice[ice < 15] = np.nan

    ice_avg.append(np.nanmean(ice))
    ice_std.append(np.nanstd(ice))
    ice_p16.append(np.nanpercentile(ice, 16))
    ice_p85.append(np.nanpercentile(ice, 86))

ice_avg = np.array(ice_avg)
ice_std = np.array(ice_std)
ice_p16 = np.array(ice_p16)
ice_p85 = np.array(ice_p85)

mpl.rcParams['xtick.labelsize'] = 14 
plt.fill_between(mms, ice_p16, ice_p85, alpha=0.5, color='gray')
plt.plot(mms, ice_avg, 'k.-')
plt.xlabel('month', fontsize=14)
plt.ylabel('average concentration, %', fontsize=14)
plt.savefig('figure_05_monthly_sic_aver.png')
plt.close()


### =====================
#####  FIGURE 7. Propagation and weighted average
# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

c = np.load('/files/sea_ice_age/osi405c_demo_archive_filled_v1/ice_drift_nh_polstere-625_multi-oi_201210011200-201210031200.nc.npz')['c']
make_map('tmp_2012-10-01', 'sic', osi_sic_dom, dst_dom,
          vmin=0, vmax=100, cmap=cm.ice, title='2012-10-01', array=c, text='A: $C_{OBS}$')
make_map('tmp_2012-10-01', '2yi', osi_sic_dom, dst_dom,
          vmin=0, vmax=100, cmap=cm.ice, array=c, text='D: $C_{2YI}$')
make_map('tmp_2012-10-01', '1yi', osi_sic_dom, dst_dom,
          vmin=100, vmax=200, cmap=cm.ice, array=c, text='G: $C_{FYI}$')

c = np.load('/files/sea_ice_age/osi405c_demo_archive_filled_v1/ice_drift_nh_polstere-625_multi-oi_201303291200-201303311200.nc.npz')['c']
make_map('tmp_2013-04-01', 'sic', osi_sic_dom, dst_dom,
          vmin=0, vmax=100, cmap=cm.ice, title='2013-04-01', array=c, text='B: $C_{OBS}$')

c = np.load('/files/sea_ice_age/osi405c_demo_archive_filled_v1/ice_drift_nh_polstere-625_multi-oi_201309151200-201309171200.nc.npz')['c']
make_map('tmp_2013-09-15', 'sic', osi_sic_dom, dst_dom,
          vmin=0, vmax=100, cmap=cm.ice, title='2013-09-15', array=c, text='C: $C_{OBS}$')

sif = np.load('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-03-29_sia.npz')['sif']
make_map('tmp_2013-04-01',
         '1yi', osi_sic_dom, dst_dom,
          vmin=0, vmax=1, cmap=cm.ice, array=sif[0], text='H: $C_{FYI}$')
make_map('tmp_2013-04-01',
         '2yi', osi_sic_dom, dst_dom,
          vmin=0, vmax=1, cmap=cm.ice, array=sif[1], text='E: $C_{2YI}$')

sif = np.load('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-09-17_sia.npz')['sif']
make_map('tmp_2013-09-15',
         '2yi', osi_sic_dom, dst_dom,
          vmin=0, vmax=1, cmap=cm.ice, array=sif[1], text='I: $C_{2YI}$')
make_map('tmp_2013-09-15',
         '3yi', osi_sic_dom, dst_dom,
          vmin=0, vmax=1, cmap=cm.ice, array=sif[2], text='F: $C_{3YI}$')

"""
!montage\
    tmp_2012-10-01_sic.png tmp_2013-04-01_sic.png tmp_2013-09-15_sic.png\
    tmp_2012-10-01_2yi.png tmp_2013-04-01_2yi.png tmp_2013-09-15_3yi.png\
    tmp_2012-10-01_1yi.png tmp_2013-04-01_1yi.png tmp_2013-09-15_2yi.png\
   -tile 3x3 -geometry +0+0 figure_07_weighted_sif.png 
"""
save_legend(cm.ice, np.linspace(0,100,20), 'Sea Ice Age Fraction, %', 'figure_07_sif_legend.png')


make_map('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2012-10-01_sia.npz',
         'sia', osi_sic_dom, dst_dom,
         vmin=0, vmax=3, cmap=cm.thermal, text='J: $SIA$')

make_map('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-03-29_sia.npz',
         'sia', osi_sic_dom, dst_dom,
         vmin=0, vmax=3, cmap=cm.thermal, text='K: $SIA$')

make_map('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-09-15_sia.npz',
         'sia', osi_sic_dom, dst_dom,
         vmin=0, vmax=3, cmap=cm.thermal, text='L: $SIA$')

"""
!montage\
    /files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2012-10-01_sia.npz_sia.png\
    /files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-03-29_sia.npz_sia.png\
    /files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2013-09-15_sia.npz_sia.png\
   -tile 3x1 -geometry +0+0 figure_07_weighted_sia.png     
"""

save_legend(cm.thermal, np.linspace(0.,3.,13.), 'Sea Ice Age, years', 'figure_07_sia_legend.png', format='%2.1f')

## FIGURE 8. COMPARISON OF SIA with different methods
"""
montage\
 /files/sea_ice_age/fowler_nsidc_2016_f02/sia/2015-12-31_sia.npz_sia.png\
 /files/sea_ice_age/fowler_osi_fv1/sia/2015-12-31_sia.npz_sia.png\
 /files/sea_ice_age/nersc_osi_fv1_2017/sia/2016-01-01_sia.npz_sia.png\
 /files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/2016-01-01_sia.npz_sia.png\
 -tile 2x2 -geometry +0+0 figure_08_sia_compar_methods.png 
"""

save_legend('jet', np.linspace(0.,5.,13.), 'Sea Ice Age, years', 'figure_08_sia_legend.png', format='%2.1f')


### FIGURE 9. Comparison of MYI
### NSIDC
# source domain
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_sia_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')
# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')
text='NSIDC'
for yy in [2012, 2013, 2014, 2015]:
    sia = get_nsidc_raw_sia('/files/nsidc0611_seaice_age_v3/iceage.grid.week.%d.52.n.v3.bin' % yy)
    myi = (sia > 1).astype(float)
    wat = (sia == 0).astype(float)
    wat[wat == 0] = np.nan
    nmap = make_map('tmp_%d-12-31' % yy, 'nsidc_myi', nsidc_sia_dom, dst_dom,
             array=myi, vmin=0, vmax=1, cmap=cm.ice, text=text, title='%d-12-31' % yy,
             water=wat)
    text=None
make_map('tmp_2016-12-31', 'nsidc_myi', nsidc_sia_dom, dst_dom,
         array=myi, vmin=-1, vmax=0, cmap=cm.ice, title='2016-12-31')
make_map('tmp_2017-03-29', 'nsidc_myi', nsidc_sia_dom, dst_dom,
         array=myi, vmin=-1, vmax=0, cmap=cm.ice, title='2017-03-29')
"""
!montage\
    tmp_2012-12-31_nsidc_myi.png\
    tmp_2013-12-31_nsidc_myi.png\
    tmp_2014-12-31_nsidc_myi.png\
    tmp_2015-12-31_nsidc_myi.png\
    tmp_2016-12-31_nsidc_myi.png\
    tmp_2017-03-29_nsidc_myi.png\
    -tile 6x1 -geometry +0+0 figure_09_myi_nsidc.png 
"""


####================== NERSC    
# source domain
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')
text = 'NERSC'
for dst_date in [
    dt.datetime(2012,12,31),
    dt.datetime(2013,12,31),
    dt.datetime(2014,12,31),
    dt.datetime(2015,12,31),
    dt.datetime(2016,12,31),
    dt.datetime(2017,3,29)]:
    ifile = '/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/%s_sia.npz' % dst_date.strftime('%Y-%m-%d')
    myi = np.load(ifile)['myi']
    sia = np.load(ifile)['sia']
    wat = (sia < 0).astype(float)
    wat[wat == 0] = np.nan
    make_map(dst_date.strftime('tmp_%Y-%m-%d'), 'nersc_myi', osi_sic_dom, dst_dom, array=myi,
             vmin=0, vmax=1, cmap=cm.ice, text=text, water=wat)
    text=None

"""
!montage\
    tmp_2012-12-31_nersc_myi.png\
    tmp_2013-12-31_nersc_myi.png\
    tmp_2014-12-31_nersc_myi.png\
    tmp_2015-12-31_nersc_myi.png\
    tmp_2016-12-31_nersc_myi.png\
    tmp_2017-03-29_nersc_myi.png\
    -tile 6x1 -geometry +0+0 figure_09_myi_nersc.png 
"""


####============== BREMEN
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
brem_dom = Domain(osi_nsr, '-te -3843750 -5343750 3756250 5856250 -tr 12500 12500')

dates = ['2013-12-31', '2014-12-31', '2015-12-31']
for date in dates:
    brem_date = parse(date)
    brem_file = '/files/sea_ice_age/bremensit/MYI-NMYI-CORRECTION-%s.nc' % brem_date.strftime('%Y%m%d')
    brem_myi = nc4.Dataset(brem_file).variables['NMYI'][:].reshape(896, 608)

    make_map('tmp_%s' % date, 'bremen_myi', brem_dom, dst_dom,
             array=brem_myi, vmin=0, vmax=100, cmap=cm.ice)

dates = ['2012-12-31', '2016-12-31', '2017-03-29']
text = 'BREMEN'
for date in dates:
    make_map('tmp_%s' % date, 'bremen_myi', brem_dom, dst_dom,
         array=np.ones(brem_myi.shape), vmin=-10, vmax=0, cmap=cm.ice, text=text)
    text = None

"""
!montage\
    tmp_2012-12-31_bremen_myi.png\
    tmp_2013-12-31_bremen_myi.png\
    tmp_2014-12-31_bremen_myi.png\
    tmp_2015-12-31_bremen_myi.png\
    tmp_2016-12-31_bremen_myi.png\
    tmp_2017-03-29_bremen_myi.png\
    -tile 6x1 -geometry +0+0 figure_09_myi_bremen.png 
"""



####================= OSISAF
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')

osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sit_dom = Domain(osi_nsr, '-te -3845000 -5345000 3755000 5855000 -tr 10000 10000')
#osi_sit_url_fmt = 'http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/type/%Y/%m/ice_type_nh_polstere-100_multi_%Y%m%d1200.nc'
tmp_dir  = './sit_osisaf_tmp/'
osi_sit_url_fmt = 'http://thredds.met.no/thredds/fileServer/osisaf/met.no/ice/type/%Y/%m/ice_type_nh_polstere-100_multi_%Y%m%d1200.nc'

text = 'OSI-SAF'
for osi_date in [
    dt.datetime(2012,12,31),
    dt.datetime(2013,12,31),
    dt.datetime(2014,12,31),
    dt.datetime(2015,12,31),
    dt.datetime(2016,12,31),
    dt.datetime(2017,3,29)]:

    url = osi_date.strftime(osi_sit_url_fmt)
    print url
    ofile = wget.detect_filename(url)
    if not os.path.exists(os.path.join(tmp_dir, ofile)):
        wget.download(url, tmp_dir)

    n = Nansat(os.path.join(tmp_dir, ofile), mapperName='generic')
    ice_type = n['ice_type'].astype(float)
    ice_type[ice_type == 255] = np.nan
    wat = (ice_type == 1).astype(float)
    wat[wat == 0] = np.nan
    make_map(osi_date.strftime('tmp_%Y-%m-%d'), 'osisaf_myi', osi_sit_dom, dst_dom,
             array=ice_type, vmin=2, vmax=3, cmap=cm.ice, text=text, water=wat)
    text = None

"""
!montage\
    tmp_2012-12-31_osisaf_myi.png\
    tmp_2013-12-31_osisaf_myi.png\
    tmp_2014-12-31_osisaf_myi.png\
    tmp_2015-12-31_osisaf_myi.png\
    tmp_2016-12-31_osisaf_myi.png\
    tmp_2017-03-29_osisaf_myi.png\
    -tile 6x1 -geometry +0+0 figure_09_myi_osisaf.png 
"""

save_legend(cm.ice, np.linspace(0,100,20), 'Multi-Year Ice Concentration, %', 'figure_09_myi_legend.png')

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
