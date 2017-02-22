import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *

from iceagelib import *

def get_mean_age(idir, thedate):
    ### weighted average of fractional ice age
    rfiles = sorted(glob.glob(idir + '*%s.npz' % thedate.strftime('%y-%m-%d')), reverse=True)
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
    
    return ice_age_weights, ice_age_mean

def add_fyi(ice_age_weights, ice_age_mean, ice_mask):
    myi_weight = ice_age_weights.sum(axis=0)
    fyi_weight = 1 - myi_weight
    fyi_weight[ice_mask == 0] = 0
    
    ice_age_weights = np.array([fyi_weight] + list(ice_age_weights))
    ice_age_weighted = np.multiply(range(1,1+len(ice_age_weights)), ice_age_weights.T).T
    ice_age_mean = ice_age_weighted.sum(axis=0) / ice_age_weights.sum(axis=0)
    
    ice_age_mean[ice_mask == 0] = 0
    ice_age_mean[np.isnan(ice_mask)] = np.nan
    
    return ice_age_weights, ice_age_mean
    
osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'
osi_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 12500 12500')

nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')
    

# get NCIDC ice age for 2015, 22 (2015-06-04)
idir_ia = '/files/nsidc0611_seaice_age_v3/'
nsidc_f = sorted(glob.glob(idir_ia + '*.bin'))[-1]
nsidc_age = np.fromfile(nsidc_f, np.uint8).reshape(361*2,361*2).astype(np.float32)
nsidc_age[nsidc_age==255] = np.nan
nsidc_age /= 5.
ice_mask = (nsidc_age > 0).astype(np.float32)
ice_mask[np.isnan(nsidc_age)] = np.nan

agew_nsidc, agem_nsidc = get_mean_age(
                            '/files/sea_ice_age/nsidc_f2_newprop_2015/',
                             dt.datetime(2015,06,4))

agew_nsidc, agem_nsidc = add_fyi(agew_nsidc, agem_nsidc, ice_mask)

agew_osi, agem_osi = get_mean_age(
                             '/files/sea_ice_age/osi_newprop_f5_zoom1/',
                             dt.datetime(2015,06,4))

agem_osi_pro = reproject_ice(osi_dom, nsidc_dom, agem_osi)
agew_osi_pro = np.zeros((2,) + agem_osi_pro.shape)
agew_osi_pro[0] = reproject_ice(osi_dom, nsidc_dom, agew_osi[0])
agew_osi_pro[1] = reproject_ice(osi_dom, nsidc_dom, agew_osi[1])

agew_osi_pro, agem_osi_pro = add_fyi(agew_osi_pro, agem_osi_pro, ice_mask)

#plt.subplot(1,2,1);plt.imshow(agem_nsidc)
#plt.subplot(1,2,2);plt.imshow(agem_osi_pro)
#plt.show()

fig = plt.figure(1, figsize=(15,5))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222, sharex=ax1, sharey=ax1)
ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1)
ax4 = fig.add_subplot(224)

def click(event):
    if event.button==1 and event.inaxes:#select initial conditions pressing left mouse button 
        x0,y0 = event.xdata,event.ydata
        print x0,y0
        aw_nsidc = agew_nsidc[:, y0, x0]
        aw_osi = agew_osi_pro[:, y0, x0]
        am_nsidc = agem_nsidc[y0, x0]
        am_osi = agem_osi_pro[y0, x0]
        ice_age_nsidc = nsidc_age[y0, x0]
        wx = range(1,1+len(aw_nsidc))
        print wx, aw_nsidc
        ax4.clear()

        ax4.plot(wx, aw_nsidc, '.-r')
        ax4.plot([am_nsidc, am_nsidc], [0, 1], 'r-')

        ax4.plot(wx, aw_osi, '.-g')
        ax4.plot([am_osi, am_osi], [0, 1], 'g-')
        
        ax4.plot([ice_age_nsidc, ice_age_nsidc], [0, 1], 'b-')

        fig.canvas.draw()

ims1 = ax1.imshow(agem_nsidc, cmap='jet', vmax=3)#;plt.colorbar(ims1)
ims2 = ax2.imshow(agem_osi_pro, cmap='jet', vmax=3)#;plt.colorbar(ims2)
ims3 = ax3.imshow(nsidc_age, cmap='jet', vmax=3)#;plt.colorbar(ims2)
plt.connect('button_press_event', click)
plt.show()

water = np.ma.array(np.zeros_like(nsidc_age), mask=nsidc_age != 0)

plt.close('all')
r1,r2 = 150,540
c1,c2 = 180,540
fig,axs = plt.subplots(1,3)
ax0 = axs[0].imshow(nsidc_age[r1:r2, c1:c2], cmap='jet', vmin=1, vmax=3)
ax0wm = axs[0].imshow(water[r1:r2, c1:c2], cmap='gray')
axs[0].set_title('Ice Age from NSIDC', fontsize=8)

ax1 = axs[1].imshow(agem_nsidc[r1:r2, c1:c2], cmap='jet', vmin=1, vmax=3)
ax1wm = axs[1].imshow(water[r1:r2, c1:c2], cmap='gray')
axs[1].set_title('NSIDC drift, NERSC ice age', fontsize=8)

ax2 = axs[2].imshow(agem_osi_pro[r1:r2, c1:c2], cmap='jet', vmin=1, vmax=3)
ax2wm = axs[2].imshow(water[r1:r2, c1:c2], cmap='gray')
axs[2].set_title('OSI drift, NERSC ice age', fontsize=8)

[(ax.set_xticks([]),ax.set_yticks([])) for ax in axs]
cbaxes = fig.add_axes([0.1, 0.25, 0.8, 0.03]) 
cb = fig.colorbar(ax2, cax = cbaxes, orientation='horizontal')
cb.ax.tick_params(labelsize=8)
cb.set_label('Sea Ice Age, 04 June 2015', size=8)
plt.tight_layout()
plt.savefig('new_age_maps.png', dpi=200, bbox_inches='tight', pad_inches=0)
plt.close()



raise
"""
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
    
"""
