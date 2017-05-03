import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *

from iceagelib import *
osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'

nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 6250 6250')

sit_dom = Domain(osi_nsr, '-te -3845000 -5345000 3755000 5855000 -tr 10000 10000')
sit_url_fmt = 'http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/type/%Y/%m/ice_type_nh_polstere-100_multi_%Y%m%d1200.nc'

odir = 'ages2015_conc_f10/'

# get NCIDC ice age
idir_ia = '/files/nsidc0611_seaice_age_v3/'
nsidc_files = sorted(glob.glob(idir_ia + 'iceage.grid.week.201[2,3,4,5].??.n.v3.bin'))

osimyi_list = []
for nsidc_f in nsidc_files:
    nsidc_date = get_nsidc_date(nsidc_f)
    print nsidc_date

    sit = Dataset(nsidc_date.strftime(sit_url_fmt)).variables['ice_type'][0]
    sit_pro = reproject_ice(sit_dom, nsidc_dom, sit).astype('float32')
    sit_pro[np.isnan(ice_mask)] = np.nan
    sit_pro[sit_pro == 1] = 0
    sit_pro[sit_pro == 2] = 0
    sit_pro[sit_pro == 3] = 1
    sit_pro[sit_pro == 4] = 0.5
    sit_pro[sit_pro == 255] = np.nan

    osisaf_sit_list.append(sit_pro)

