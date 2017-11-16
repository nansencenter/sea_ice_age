import json

import numpy as np
import matplotlib.pyplot as plt
import ogr

from nansat import *
from iceagelib import *
import netCDF4 as nc4


# read SAR mosaic
factor = 4
m = plt.imread('/mnt/data/nersc/papers/1702_ice_age/primo2016SAR.png')[::factor, ::factor]
print m.shape

## READ Leif outline
shpFile = ogr.Open('MY_main_20160101.shp')
shape = shpFile.GetLayer(0)
feature_count = shape.GetFeatureCount()
coords = [json.loads(shape.GetFeature(i).ExportToJson())['geometry']['coordinates']
          for i in range(feature_count)]
shpFile = None

# mosaic domain
res = 1000
mos_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=-45  +lat_ts=70 +no_defs')
#mos_dom = Domain(mos_nsr, '-te -2350000 -3400000 2150000 2100000 -tr %d %d' % (res*factor, res*factor))
width = m.shape[1] * res * factor
height = m.shape[0] * res * factor
minx = -2350000
miny = -3400000
maxx = minx + width
maxy = miny + height
mos_dom = Domain(mos_nsr, '-te %d %d %d %d -tr %d %d' % (minx, miny, maxx, maxy, res*factor, res*factor))
print mos_dom.shape()

# NSIDC
'''
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_sia_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')
nsidc_file = '/files/nsidc0611_seaice_age_v3/iceage.grid.week.2015.52.n.v3.bin'
nsidc_sia = get_nsidc_raw_sia(nsidc_file)
nsidc_myi = gaussian_filter((nsidc_sia > 1).astype(float), 2)
nsidc_myi_pro = reproject_ice(nsidc_sia_dom, mos_dom, nsidc_myi)
'''

# SICCI
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')
ifile = '/files/sea_ice_age/nersc_osi_fv6_2017_conc/sia/2016-01-01_sia.npz'
sicci_myi = np.load(ifile)['myi']
sicci_myi_pro = reproject_ice(osi_sic_dom, mos_dom, sicci_myi)

# BREMEN
'''
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
breme_dom = Domain(osi_nsr, '-te -3843750 -5343750 3756250 5856250 -tr 12500 12500')
breme_file = '/files/sea_ice_age/bremensit/MYI-NMYI-CORRECTION-20160101.nc'
breme_myi = nc4.Dataset(breme_file).variables['NMYI'][:].reshape(896, 608) / 100.
breme_myi[np.isnan(breme_myi)] = 0
breme_myi = gaussian_filter(breme_myi, 1)
breme_myi_pro = reproject_ice(breme_dom, mos_dom, breme_myi)
'''


plt.close('all')
nmap = Nansatmap(mos_dom)
nmap.imshow(m, vmin=0.15, vmax=0.85, cmap='gray')
for seg_coords in coords:
    seg_coords = np.array(seg_coords)*1.03 + np.array([80000, 90000])
    nmap.plot(seg_coords.T[0]-minx, seg_coords.T[1]-miny, '-', color='#00ffff', linewidth=0.3)

levels=[0.1, 0.5]
linewidths=[0.7, 0.3]
nmap.contour(sicci_myi_pro, levels=levels, colors='r', linewidths=linewidths, label=False)
#plt.contour(breme_myi_pro, levels=levels, colors='y', linewidths=linewidths)
#plt.contour(nsidc_myi_pro, levels=levels, colors='b', linewidths=linewidths)
#for seg_coords in coords:
#    seg_coords = np.array(seg_coords) * 1.03 + np.array([100000, 70000])
#    nmap.plot(seg_coords.T[0]+2350000, seg_coords.T[1]+3400000, '-b', linewidth=0.1)


#levels=[0.8, 1]
#alpha=0.1
#plt.contourf(sicci_myi_pro, levels=levels, colors='r', linewidths=0, alpha=alpha)
#plt.contourf(nsidc_myi_pro, levels=levels, colors='b', linewidths=0, alpha=alpha)

#plt.show()
plt.xticks([])
plt.yticks([])
plt.ylim([1500000, 5300000])
plt.xlim([300000, 4000000])
plt.savefig('figure_12_SAR_mosaic_fv6.png', pad_inches=0, bbox_inches='tight', dpi=600)
plt.close()
