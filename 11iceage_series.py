import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *

from iceagelib import *

osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'
osi_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 12500 12500')

nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')

sit_dom = Domain(osi_nsr, '-te -3845000 -5345000 3755000 5855000 -tr 10000 10000')
sit_url_fmt = 'http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/type/%Y/%m/ice_type_nh_polstere-100_multi_%Y%m%d1200.nc'

odir = 'ages2015_conc/'

# get NCIDC ice age for 2015, 22 (2015-06-04)
idir_ia = '/files/nsidc0611_seaice_age_v3/'
nsidc_files = sorted(glob.glob(idir_ia + 'iceage.grid.week.2015.??.n.v3.bin'))
k = 0
for nsidc_f in nsidc_files:
    nsidc_date = get_nsidc_date(nsidc_f)
    print nsidc_date

    nsidc_age = np.fromfile(nsidc_f, np.uint8).reshape(361*2,361*2).astype(np.float32)
    nsidc_age[nsidc_age==255] = np.nan
    nsidc_age /= 5.
    ice_mask = (nsidc_age > 0).astype(np.float32)
    ice_mask[np.isnan(nsidc_age)] = np.nan

    #agew_nsidc, agem_nsidc = get_mean_age(
    #                            '/files/sea_ice_age/nsidc_f2_newprop_2015/',
    #                             nsidc_date)
    #agew_nsidc, agem_nsidc = add_fyi(agew_nsidc, agem_nsidc, ice_mask)

    agew_osi, agem_osi = get_mean_age(
                                 '/files/sea_ice_age/osi_newprop_f5_zoom1_conc/',
                                 nsidc_date)

    agem_osi_pro = reproject_ice(osi_dom, nsidc_dom, agem_osi)
    agew_osi_pro = np.zeros((agew_osi.shape[0],) + agem_osi_pro.shape)
    for i in range(agew_osi.shape[0]):
        agew_osi_pro[i] = reproject_ice(osi_dom, nsidc_dom, agew_osi[i])

    agew_osi_pro, agem_osi_pro = add_fyi(agew_osi_pro, agem_osi_pro, ice_mask)

    # read sea ice type from OSI-SAF and reproject
    sit = Dataset(nsidc_date.strftime(sit_url_fmt)).variables['ice_type'][0]
    sit_pro = reproject_ice(sit_dom, nsidc_dom, sit).astype('float32')
    sit_pro[np.isnan(ice_mask)] = np.nan
    
    water = np.ma.array(np.zeros_like(nsidc_age), mask=nsidc_age != 0)

    plt.close('all')
    r1,r2 = 150,540
    c1,c2 = 180,540
    fig,axs = plt.subplots(2, 3)

    ax0 = axs[0,0].imshow(nsidc_age[r1:r2, c1:c2], cmap='jet', vmin=1, vmax=4)
    ax0wm = axs[0,0].imshow(water[r1:r2, c1:c2], cmap='gray')
    axs[0,0].set_title('SIA NSIDC', fontsize=8)

    #ax1 = axs[1].imshow(agem_nsidc[r1:r2, c1:c2], cmap='jet', vmin=1, vmax=3)
    #ax1wm = axs[1].imshow(water[r1:r2, c1:c2], cmap='gray')
    #axs[1].set_title('NSIDC drift, NERSC ice age', fontsize=8)

    ax2 = axs[0,1].imshow(agem_osi_pro[r1:r2, c1:c2], cmap='jet', vmin=1, vmax=4)
    ax2wm = axs[0,1].imshow(water[r1:r2, c1:c2], cmap='gray')
    axs[0,1].set_title('SIA NERSC', fontsize=8)

    ax1 = axs[0,2].imshow(sit_pro[r1:r2, c1:c2], cmap='jet', vmin=1, vmax=4)
    ax1wm = axs[0,2].imshow(water[r1:r2, c1:c2], cmap='gray')
    axs[0,2].set_title('SIT OSISAF', fontsize=8)

    age2 = agew_osi_pro[1]
    age2[np.isnan(ice_mask)] = np.nan
    ax3 = axs[1,0].imshow(age2[r1:r2, c1:c2], cmap='jet', vmin=0, vmax=1)
    ax3wm = axs[1,0].imshow(water[r1:r2, c1:c2], cmap='gray')
    axs[1,0].set_title('2YO NERSC', fontsize=8)

    age3 = agew_osi_pro[2]
    age3[np.isnan(ice_mask)] = np.nan
    ax3 = axs[1,1].imshow(age3[r1:r2, c1:c2], cmap='jet', vmin=0, vmax=1)
    ax3wm = axs[1,1].imshow(water[r1:r2, c1:c2], cmap='gray')
    axs[1,1].set_title('3YO NERSC', fontsize=8)

    age4 = agew_osi_pro[3]
    age4[np.isnan(ice_mask)] = np.nan
    ax3 = axs[1,2].imshow(age4[r1:r2, c1:c2], cmap='jet', vmin=0, vmax=1)
    ax3wm = axs[1,2].imshow(water[r1:r2, c1:c2], cmap='gray')
    axs[1,2].set_title('4YO NERSC', fontsize=8)

    [(ax.set_xticks([]),ax.set_yticks([])) for ax in axs.flatten()]

    cbaxes = fig.add_axes([0.0, 0, 0.45, 0.03]) 
    cb = fig.colorbar(ax2, cax=cbaxes, orientation='horizontal')
    cb.ax.tick_params(labelsize=8)
    cb.set_label('SIA, %s' % nsidc_date.strftime('%Y-%m-%d'), size=8)

    cbaxes2 = fig.add_axes([0.5, 0, 0.45, 0.03]) 
    cb2 = fig.colorbar(ax3, cax=cbaxes2, orientation='horizontal')
    cb2.ax.tick_params(labelsize=8)
    cb2.set_label('MYI fraction, %s' % nsidc_date.strftime('%Y-%m-%d'), size=8)

    plt.tight_layout()
    plt.savefig(odir + 'new_age_maps_%03d_conc.png' % k, dpi=200, bbox_inches='tight', pad_inches=0)
    plt.close()
    k += 1
