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


### FIGURE 9. Comparison of MYI
# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -2000000 -2400000 2000000 2600000 -tr 10000 10000')


### NSIDC
# source domain
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_sia_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')
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
text = 'SICCI'
for dst_date in [
    dt.datetime(2012,12,31),
    dt.datetime(2013,12,31),
    dt.datetime(2014,12,31),
    dt.datetime(2015,12,31),
    dt.datetime(2016,12,31),
    dt.datetime(2017,3,29)]:
    ifile = '/files/sea_ice_age/nersc_osi_fv2_2017_conc_repro/sia/%s_sia.npz' % dst_date.strftime('%Y-%m-%d')
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

dates = ['2013-12-31', '2014-12-31', '2015-12-31', '2016-12-31', '2017-03-29']
for date in dates:
    brem_date = parse(date)
    brem_file = '/files/sea_ice_age/bremensit/MYI-NMYI-CORRECTION-%s.nc' % brem_date.strftime('%Y%m%d')
    brem_myi = nc4.Dataset(brem_file).variables['NMYI'][:].reshape(896, 608)

    make_map('tmp_%s' % date, 'bremen_myi', brem_dom, dst_dom,
             array=brem_myi, vmin=0, vmax=100, cmap=cm.ice)

dates = ['2012-12-31']
text = 'BU'
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
!montage\
    figure_09_myi_nsidc.png\
    figure_09_myi_nersc.png\
    figure_09_myi_osisaf.png\
    figure_09_myi_bremen.png\
    -tile 1x4 -geometry +0+0 figure_09_myi.png 
"""
