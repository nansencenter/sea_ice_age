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

import pandas as pd


#### SHOW FOR 2015
"""
# destination domain
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
dst_dom = Domain(dst_nsr, '-te -900000 -900000 900000 900000 -tr 10000 10000')

#ifiles = sorted(glob.glob('iabp/*/*.csv')) + sorted(glob.glob('iabp/*/*.dat'))
ifiles = sorted(glob.glob('iabp/201[2,3,4,5]/*.dat'))

nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_sia_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')

bad_files = [
'iabp/2015/300234061261790.dat',
'iabp/2015/300234062429080.dat',
'iabp/2015/300234061263820.dat',
'iabp/2015/300234062427080.dat',
'iabp/2015/300234061264830.dat',
'iabp/2015/300234061266810.dat',
]

for year in [2015]:
    for wn in [36, 37, 38, 39, 40, 41]:
        print year, wn
        doy_min = wn * 7 - 15
        doy_max = wn * 7

        sia = get_nsidc_raw_sia('/files/nsidc0611_seaice_age_v3/iceage.grid.week.%04d.%02d.n.v3.bin' % (year, wn))
        sia_pro = reproject_ice(nsidc_sia_dom, dst_dom, sia).astype('float32')
        nmap = Nansatmap(dst_dom)
        nmap.imshow(sia_pro, cmap='plasma', vmax=5)
        lons = []
        lats = []
        ifls = []
        for ifile in ifiles:
            if ifile in bad_files:
                continue
            names = open(ifile).readline().strip().split(',')
            if ifile.endswith('.dat'):
                data = pd.read_table(ifile, skiprows=1, names=names, delim_whitespace=True)
            else:
                data = pd.read_table(ifile, skiprows=1, names=names, sep=',')
            lon, lat, years, doy = map(np.array, [data['Lon'], data['Lat'], data['Year'], data['DOY']])
            gpi = np.isfinite(lon) * np.isfinite(lat) * (years == year) * (doy >= doy_min) * (doy <= doy_max)
            lon = lon[gpi]
            lat = lat[gpi]
            if len(lon) == 0:
                continue
            if lat.max() < 85:
                continue
            lons.append(lon)
            lats.append(lat)
            ifls.append(ifile)
            x, y = nmap(lon, lat)
            nmap.plot(x, y, '.-', ms=4)
            nmap.plot(x[-1], y[-1], '.k', ms=2)
            #plt.text(x[-1], y[-1], os.path.basename(ifile))
        nmap.draw_continents()
        nmap.save('sia_buoys_%04d.%02d.png' % (year, wn), dpi=300)
        plt.close('all')
#"""

#### SHOW for 1979

ifile = '/files/sea_ice_age/fowler_nsidc_1985_f04/icemap_1979-09-03_1979-10-22.npz'
ice = np.load(ifile)['ice']
ice += 1
ice[np.isnan(ice)] = 0
# source domain
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_sia_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 6250 6250')
nsidc_ice = Nansat(domain=nsidc_sia_dom, array=ice)

nsidc_ice.crop(461,420, 500,500)
ice = nsidc_ice[1]

nmap = Nansatmap(nsidc_ice, resolution='l', lat_0=90, lon_0=0)
nmap.imshow(ice, cmap=cm.ice, vmin=0, vmax=2)
nmap.draw_continents()

year = 1979
doy_min = 264
doy_max = 294

ifiles = sorted(glob.glob('iabp/1979/*.csv'))
for ifile in ifiles:
    names = open(ifile).readline().strip().split(',')
    data = pd.read_table(ifile, skiprows=1, names=names, sep=',')
    lon, lat, years, doy = map(np.array, [data['Lon'], data['Lat'], data['Year'], data['DOY']])
    gpi = np.isfinite(lon) * np.isfinite(lat) * (years == year) * (doy >= doy_min) * (doy <= doy_max)
    lon = lon[gpi]
    lat = lat[gpi]
    if len(lon) == 0:
        continue
    if lat.max() < 85:
        continue
    x, y = nmap(lon, lat)
    nmap.plot(x, y, '.-r', ms=4)
    nmap.plot(x[-1], y[-1], '.k', ms=2)
linewidth=0.5
nmap.drawmeridians(range(-180, 180, 30), (1,1,1,1), linewidth=linewidth, labels=[False, False, False, True, ])
nmap.drawparallels(range(70, 90, 5), (1,1,1,1), linewidth=linewidth, )
nmap.save('sia_buoys_1979-10-22.png', dpi=300)
plt.close('all')
