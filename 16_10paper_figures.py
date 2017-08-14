import os
import glob
import datetime as dt
from dateutil.parser import parse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

import wget
import netCDF4 as nc4
import wget

from ovl_plugins.lib.interpolation import fill_gaps_nn
from nansat import *
from iceagelib import *

from cmocean import cm

### FIGURE 10: Time series
# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')

age_range = range(10)

### NSIDC
"""
# source domain
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_sia_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')

nsidc_files = sorted(glob.glob('/files/nsidc0611_seaice_age_v3/iceage.grid.week.201[2,3,4,5,6].*.n.v3.bin'))
#nsidc_files = sorted(glob.glob('/files/nsidc0611_seaice_age_v3/iceage.grid.week.2012.*.n.v3.bin'))
nsidc_dates = map(get_nsidc_date, nsidc_files)

nsidc_age_list = []
for nsidc_file in nsidc_files:
    print os.path.basename(nsidc_file)
    sia = reproject_ice(nsidc_sia_dom, dst_dom, get_nsidc_raw_sia(nsidc_file))
    nsidc_age_list.append(sia)

### USE the MOST CRUDE LANDMASK from NSIDC
landmask = np.isnan(sia)
np.savez('landmask', landmask=landmask)


nsidc_age_list = np.array(nsidc_age_list)
nsidc_age_list[np.isnan(nsidc_age_list)] = -1

nsidc_age_area = []
for age in age_range:
    nsidc_age_area.append((nsidc_age_list == age).sum(axis=(1,2)))
nsidc_age_area = np.array(nsidc_age_area)
np.savez('nsidc_age_area', nsidc_age_area=nsidc_age_area, nsidc_dates=nsidc_dates)
"""

# NERSC
"""
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')
landmask = np.load('landmask.npz')['landmask']

nersc_sia_files = sorted(glob.glob('/files/sea_ice_age/nersc_osi_fv1_2017_conc/sia/201[2,3,4,5,6,7]*_sia.npz'))
nersc_age_area = []
nersc_dates = []
for nersc_sia_file in nersc_sia_files:
    print os.path.basename(nersc_sia_file)
    nersc_dates.append(parse(os.path.basename(nersc_sia_file).split('_')[0]))
    sif = np.load(nersc_sia_file)['sif']
    sif_pro = np.array([reproject_ice(osi_sic_dom, dst_dom, sif_n) for sif_n in sif])
    # replace vallid pixels over land with Nan (use crude landmask)
    sif_pro[:, landmask] = np.nan
    # replace NaN pixels over water (crude land mask) with 0
    sif_pro[:, np.isnan(sif_pro[0]) * (~landmask)] = 0

    noland_area = np.isfinite(sif_pro[0]).sum()
    ice_area = np.nansum(sif_pro)
    water_area = noland_area - ice_area

    age_area = np.zeros(len(age_range))
    age_area[0] = water_area
    for age in age_range:
        if age >= sif_pro.shape[0]:
            break
        age_area[age+1] = np.nansum(sif_pro[age])
    
    nersc_age_area.append(age_area)
nersc_age_area = np.array(nersc_age_area).T
np.savez('nersc_age_area', nersc_age_area=nersc_age_area, nersc_dates=nersc_dates)
#"""


### BREMEN
"""
landmask = np.load('landmask.npz')['landmask']
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
brem_dom = Domain(osi_nsr, '-te -3843750 -5343750 3756250 5856250 -tr 12500 12500')
brem_files = sorted(glob.glob('/files/sea_ice_age/bremensit/MYI-NMYI-CORRECTION-*.nc'))
brem_myi_areas = []
brem_dates = []
for brem_file in brem_files:
    print os.path.basename(brem_file)
    brem_dates.append(parse(os.path.splitext(brem_file)[0].split('-')[3]))
    brem_myi = nc4.Dataset(brem_file).variables['NMYI'][:].reshape(896, 608)
    brem_myi_pro = reproject_ice(brem_dom, dst_dom, brem_myi)
    # replace vallid pixels over land with Nan (use crude landmask)
    brem_myi_pro[landmask] = np.nan
    # replace NaN pixels over water (crude land mask) with 0
    brem_myi_pro[np.isnan(brem_myi_pro) * (~landmask)] = 0

    brem_myi_areas.append(np.nansum(brem_myi_pro/100.))
brem_myi_areas = np.array(brem_myi_areas)
np.savez('brem_myi_areas.npz', brem_myi_areas=brem_myi_areas, brem_dates=brem_dates)
#"""


### OSI SAF SIT
"""
landmask = np.load('landmask.npz')['landmask']
osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'
sit_dom = Domain(osi_nsr, '-te -3845000 -5345000 3755000 5855000 -tr 10000 10000')
sit_url_fmt = 'http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/type/%Y/%m/ice_type_nh_polstere-100_multi_%Y%m%d1200.nc'
sit_url_fmt = '/files/sea_ice_age/osi_ice_type_nh/ice_type_nh_polstere-100_multi_%Y%m%d1200.nc'

osi_dates_tmp = pd.date_range('2012-01-01', '2016-12-31', freq='1D')
osi_myi_areas = []
osi_dates = []
for osi_date in osi_dates_tmp:
    print osi_date
    fname = osi_date.strftime(sit_url_fmt)
    if not os.path.exists(fname):
        continue
    ds = nc4.Dataset(fname)
    sit = ds.variables['ice_type'][0]
    flg = ds.variables['status_flag'][0]
    msk = (flg == 101) + (flg == 102)
    sitf = fill_gaps_nn(sit, 50, mask=msk)
    sit_pro = reproject_ice(sit_dom, dst_dom, sitf).astype('float32')
    sit_pro[sit_pro != 3] = 0
    sit_pro[sit_pro == 3] = 1
    sit_pro[landmask] = np.nan

    osi_myi_areas.append(np.nansum(sit_pro))
    osi_dates.append(osi_date)

osi_doys = np.array([(d0- dt.datetime(d0.year,1,1)).days for d0 in osi_dates])
osi_gpi = (osi_doys > 258) + (osi_doys < 135)
np.savez('osi_myi_areas.npz', osi_myi_areas=osi_myi_areas, osi_dates=osi_dates, osi_gpi=osi_gpi)
#"""


### PLOT
nsidc_age_area = np.load('nsidc_age_area.npz')['nsidc_age_area']
nsidc_dates = np.load('nsidc_age_area.npz')['nsidc_dates']
nsidc_age_cum_area = np.cumsum(nsidc_age_area[::-1], axis=0)[::-1]

nersc_age_area = np.load('nersc_age_area.npz')['nersc_age_area']
nersc_dates = np.load('nersc_age_area.npz')['nersc_dates']
nersc_age_cum_area = np.cumsum(nersc_age_area[::-1], axis=0)[::-1]

osi_myi_areas = np.load('osi_myi_areas.npz')['osi_myi_areas']
osi_dates = np.load('osi_myi_areas.npz')['osi_dates']
osi_gpi = np.load('osi_myi_areas.npz')['osi_gpi']

brem_myi_areas = np.load('brem_myi_areas.npz')['brem_myi_areas']
brem_dates = np.load('brem_myi_areas.npz')['brem_dates']

ax1 = plt.subplot(2,1,1)
for i in range(1, nsidc_age_cum_area.shape[0]):
    plt.fill_between(nsidc_dates, nsidc_age_cum_area[i]*100/1000000., 0, label='%s-y.o.' % (i-1))
plt.plot(brem_dates, brem_myi_areas*100/1000000., 'k.', ms=3)
plt.plot(osi_dates[osi_gpi], osi_myi_areas[osi_gpi]*100/1000000., '*', ms=1, c='#ffff00')

#plt.title('NSIDC')
plt.ylabel('Area, 10$^6$ km$^2$')
plt.ylim([0,8.5])
plt.setp(ax1.get_xticklabels(), visible=False)
#plt.legend(['Bremen', 'OSI', '0-YO', '1-YO', '2-YO', '3-YO', '4-YO', '5-YO', ], loc=2, fontsize=8)
plt.legend(['UoB', 'OSI-SAF', '1st YI', '2nd YI', '3rd YI', '4th YI', '5th YI', '6th YI', ],
           loc=2, fontsize=8, ncol=2)
plt.text(dt.datetime(2015,8,1), 7.5, 'NSIDC', fontsize=14)

ax2 = plt.subplot(2,1,2, sharex=ax1)
for i in range(1, nersc_age_cum_area.shape[0]):
    plt.fill_between(nersc_dates, nersc_age_cum_area[i]*100/1000000., 0, label='%s-y.o.' % (i-1))
plt.plot(brem_dates, brem_myi_areas*100/1000000., 'k.', ms=3)
plt.plot(osi_dates[osi_gpi], osi_myi_areas[osi_gpi]*100/1000000., '*', ms=1, c='#ffff00')
plt.text(dt.datetime(2015,8,1), 7.5, 'NERSC', fontsize=14)

#plt.title('NERSC')
plt.ylabel('Area, 10$^6$ km$^2$')
plt.ylim([0,8.5])
#plt.xlim([dt.datetime(2013,1,1), dt.datetime(2016,1,1)])
plt.xlabel('Date')

plt.tight_layout(pad=0)
plt.savefig('figure_10_sia_components_ts.png', dpi=300, )
plt.close('all')
