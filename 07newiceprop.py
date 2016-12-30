### RUN ADVECTION OF ICE AGE (FROM NSIDC) USING OVL_PLUGINS

import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from ovl_plugins.lib.lagrangian import rungekutta4
from ovl_plugins.propagate.propagate import propagate_product_xy

from iceagelib import *

idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
idir_ia = '/files/nsidc0611_seaice_age_v3/'
odir = '/files/sea_ice_age/nsidc_f2_newprop/'
res = 12500
factor = 2

ifiles = sorted(glob.glob(idir_uv + 'icemotion.grid.week.19*.bin'))

"""
propagate_from_newprop(get_i_of_file_nsidc(1978, 44, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop/icemotion.grid.week.1978.44')

propagate_from_newprop(get_i_of_file_nsidc(1979, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop/icemotion.grid.week.1979.35')

propagate_from_newprop(get_i_of_file_nsidc(1980, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop/icemotion.grid.week.1980.35')

propagate_from_newprop(get_i_of_file_nsidc(1981, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop/icemotion.grid.week.1981.35')

propagate_from_newprop(get_i_of_file_nsidc(1982, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop/icemotion.grid.week.1982.35')

propagate_from_newprop(get_i_of_file_nsidc(1983, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop/icemotion.grid.week.1983.35')

propagate_from_newprop(get_i_of_file_nsidc(1984, 35, ifiles), ifiles, res=res, factor=factor, odir=odir)
vis_ice_npz('/files/sea_ice_age/nsidc_f2_newprop/icemotion.grid.week.1984.35')
"""

# get NCIDC ice age for 1984, 52
nsidc_f = sorted(glob.glob(idir_ia + '*.bin'))[51]
nsidc_age = np.fromfile(nsidc_f, np.uint8).reshape(361*2,361*2).astype(np.float32)
nsidc_age[nsidc_age==255] = np.nan
nsidc_age /= 5.

### weighted average of fractional ice age
rfiles = sorted(glob.glob(odir + '*%04d.%02d.n.v3.bin_icemap.npy.npz' % (1984, 52)), reverse=True)
age0 = np.load(rfiles[0])['ice']
ice_age_sum = np.zeros(age0.shape)
ice_age_wsum = np.zeros(age0.shape)
ice_age_weights = []
theage = 2
for rfile in rfiles:
    weight = np.load(rfile)['ice']
    ice_age_weights.append(weight)
    age = np.zeros(weight.shape) + theage
    ice_age_sum += age * weight
    ice_age_wsum += weight
    theage += 1

ice_age_weights = np.array(ice_age_weights)

ice_age_mean = ice_age_sum / ice_age_wsum
# add one year ice and water/landmask from nsidc_age
ice_age_mean[np.isnan(ice_age_mean)] = 1
#ice_age_mean[nsidc_age == 0] = 0
ice_age_mean[np.isnan(nsidc_age)] = np.nan

fig = plt.figure(1, figsize=(15,5))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132, sharex=ax1, sharey=ax1)
ax3 = fig.add_subplot(133)

def click(event):
    if event.button==1 and event.inaxes:#select initial conditions pressing left mouse button 
        x0,y0 = event.xdata,event.ydata
        weights = ice_age_weights[:, y0, x0]
        ice_age_nsidc = nsidc_age[y0, x0]
        mean_age = ice_age_mean[y0, x0]
        print weights
        ax3.clear()
        ax3.plot(range(2, 9), weights, '.-')
        ax3.plot([ice_age_nsidc, ice_age_nsidc], [0, 1], 'r-')
        ax3.plot([mean_age, mean_age], [0, 1], 'g-')
        fig.canvas.draw()

ims1 = ax1.imshow(nsidc_age, cmap='jet', vmax=8)#;plt.colorbar(ims1)
ims2 = ax2.imshow(ice_age_mean, cmap='jet', vmax=8)#;plt.colorbar(ims2)
plt.connect('button_press_event', click)
plt.show()

# plot histograms with areas
nsidc_hist = []
ages = np.arange(2, 9)
for age in ages:
    nsidc_hist.append(len(nsidc_age[nsidc_age == age]))
plt.bar(ages, nsidc_hist, 0.5, label='NSIDC')

nersc_hist = []
for w in ice_age_weights:
    nersc_hist.append(np.nansum(w))
plt.bar(ages+0.5, nersc_hist, 0.5, label='NERSC')
plt.legend(loc=0)
plt.show()
    

