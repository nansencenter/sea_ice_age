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

### FIGURE 11: Interannual dynamics

dates = [
dt.datetime(2013,1,1),
dt.datetime(2014,1,1),
dt.datetime(2015,1,1),
dt.datetime(2016,1,1),
]

# load data
nsidc_age_area = np.load('nsidc_age_area.npz')['nsidc_age_area']
nsidc_dates = np.load('nsidc_age_area.npz')['nsidc_dates']

nersc_age_area = np.load('nersc_age_area.npz')['nersc_age_area']
nersc_dates = np.load('nersc_age_area.npz')['nersc_dates']

osi_myi_areas = np.load('osi_myi_areas.npz')['osi_myi_areas']
osi_dates = np.load('osi_myi_areas.npz')['osi_dates']
osi_gpi = np.load('osi_myi_areas.npz')['osi_gpi']

brem_myi_areas = np.load('brem_myi_areas.npz')['brem_myi_areas']
brem_dates = np.load('brem_myi_areas.npz')['brem_dates']

# smooth MYI
brem_myi_areas_gf = gaussian_filter(brem_myi_areas,5)
osi_myi_areas_gf = gaussian_filter(osi_myi_areas,5)

# find corresponding data
nsidc_areas = []
nersc_areas = []
brm_areas = []
osi_areas = []
for date in dates:
    nsidc_areas.append(nsidc_age_area[:, np.argmin(np.abs(nsidc_dates - date))])
    nersc_areas.append(nersc_age_area[:, np.argmin(np.abs(nersc_dates - date))])
    brm_areas.append(brem_myi_areas_gf[np.argmin(np.abs(brem_dates - date))])
    osi_areas.append(osi_myi_areas_gf[np.argmin(np.abs(osi_dates - date))])

nsidc_areas = np.vstack([[0,0,0,0], np.array(nsidc_areas).T[::-1]])*100/1000000.
nersc_areas = np.vstack([[0,0,0,0], np.array(nersc_areas).T[::-1]])*100/1000000.
nsidc_areas_cs = np.cumsum(nsidc_areas, axis=0)
nersc_areas_cs = np.cumsum(nersc_areas, axis=0)

brm_areas = np.array(brm_areas)*100/1000000.
osi_areas = np.array(osi_areas)*100/1000000.

## calculate MYI reduction speed
nsidc_myi_decrease_rel = (nsidc_areas_cs[-3][:-1] - nsidc_areas_cs[-4][1:]) / nsidc_areas_cs[-3][:-1]
nsidc_myi_decrease_abs = (nsidc_areas_cs[-3][:-1] - nsidc_areas_cs[-4][1:]) 

nersc_myi_decrease_abs = nersc_areas_cs[-3][:-1] - nersc_areas_cs[-4][1:]
nersc_myi_decrease_rel = (nersc_areas_cs[-3][:-1] - nersc_areas_cs[-4][1:]) / nersc_areas_cs[-3][:-1]

print (nersc_myi_decrease_rel / nsidc_myi_decrease_rel).mean()
print (nersc_myi_decrease_abs / nsidc_myi_decrease_abs).mean()

plt.figure(figsize=(8,4))
for i in range(8,0,-1):
    plt.bar([0, 3, 6, 9], nsidc_areas[i+1], 0.95, bottom=nsidc_areas_cs[i], color='C%d' % (8-i))

plt.plot([3.5, 6.5, 9.5], brm_areas[1:], '.k')
plt.plot([0.5, 3.5, 6.5, 9.5], osi_areas, '*k', ms=7)

plt.legend(['UB', 'OSI-SAF', '1st YI', '2nd YI', '3rd YI', '4th YI', '5th YI', '6th YI', ],
           loc=1, fontsize=8)

for i in range(8,0,-1):
    plt.bar([1, 4, 7,10], nersc_areas[i+1], 0.95, bottom=nersc_areas_cs[i], color='C%d' % (8-i))

for x in [0,3,6, 9]:
    plt.text(x, 5, 'NSIDC', rotation=90)
    plt.text(x+1, 5, 'SICCI', rotation=90)

plt.xlim([-1, 12])    
plt.xticks([0.5, 3.5, 6.5, 9.5], ['2013-01-01', '2014-01-01', '2015-01-01', '2016-01-01'])
plt.xlabel('Years')
plt.ylabel('Area, 10$^6$ km$^2$')
plt.tight_layout(pad=0)
plt.savefig('figure_11_sia_interannual_fv5.png', dpi=300, )
plt.close('all')

