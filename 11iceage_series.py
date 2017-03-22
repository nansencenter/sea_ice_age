import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from nansat import *

from iceagelib import *
osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'

brem_idir = '/media/antonk/Data/files/sea_ice_age/bremensit/'
brem_dom = Domain(osi_nsr, '-te -3843750 -5343750 3756250 5856250 -tr 12500 12500')

# factor 5
osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'
osi_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 12500 12500')
osi_idir = '/files/sea_ice_age/osi_newprop_f5_zoom1_conc/'

# factor 10 # for OSI-SAF drift
factor = 4 # for NSIDC 
osi_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 6250 6250')
osi_idir = '/files/sea_ice_age/osi_newprop_f10_zoom1_conc_sigma2/'
idir_icemap_nsidc = '/files/sea_ice_age/nsidc_f4_newprop_2015/'
idir_icemap_nersc = '/files/sea_ice_age/osi_newprop_f10_zoom1_conc/'



nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 6250 6250')

sit_dom = Domain(osi_nsr, '-te -3845000 -5345000 3755000 5855000 -tr 10000 10000')
sit_url_fmt = 'http://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/type/%Y/%m/ice_type_nh_polstere-100_multi_%Y%m%d1200.nc'

odir = 'ages2015_conc_f10/'

# get NCIDC ice age for 2015, 22 (2015-06-04)
idir_ia = '/files/nsidc0611_seaice_age_v3/'
nsidc_files = sorted(glob.glob(idir_ia + 'iceage.grid.week.201[2,3,4,5].??.n.v3.bin'))

k = 0

nsidc_age_list = []
nsidc_agew_list = []
nersc_agew_list = []
osisaf_sit_list = []

doplot = True
doosisaf = False

minage = 1
maxage = 4
minc = 0
maxc = 1

for nsidc_f in nsidc_files:
    nsidc_date = get_nsidc_date(nsidc_f)
    print nsidc_date
    #import ipdb; ipdb.set_trace()

    nsidc_age = np.fromfile(nsidc_f, np.uint8).reshape(361*2,361*2).astype(np.float32)
    nsidc_age[nsidc_age==255] = np.nan
    nsidc_age /= 5.
    ice_mask = (nsidc_age > 0).astype(np.float32)
    ice_mask[np.isnan(nsidc_age)] = np.nan

    if factor > 2:
        ice_mask = zoom_nan(ice_mask, factor / 2., order=0)
        nsidc_age = zoom_nan(nsidc_age, factor / 2., order=0)

    agew_nsidc, agem_nsidc = get_mean_age(idir_icemap_nsidc, nsidc_date)
    agew_nsidc, agem_nsidc = add_fyi(agew_nsidc, agem_nsidc, ice_mask)

    try:
        agew_osi, agem_osi = get_mean_age(idir_icemap_nersc, nsidc_date)
    except:
        print 'NO OSI data!'
        continue
    
    agem_osi_pro = reproject_ice(osi_dom, nsidc_dom, agem_osi)
    agew_osi_pro = np.zeros((agew_osi.shape[0],) + agem_osi_pro.shape)
    for i in range(agew_osi.shape[0]):
        agew_osi_pro[i] = reproject_ice(osi_dom, nsidc_dom, agew_osi[i])
    agew_osi_pro, agem_osi_pro = add_fyi(agew_osi_pro, agem_osi_pro, ice_mask)
    agem_osi_pro[ice_mask != 1] = 0
    agew_osi_pro[:, ice_mask != 1] = 0
    
    # read sea ice type from OSI-SAF and reproject
    if doosisaf:
        sit = Dataset(nsidc_date.strftime(sit_url_fmt)).variables['ice_type'][0]
        sit_pro = reproject_ice(sit_dom, nsidc_dom, sit).astype('float32')
        sit_pro[np.isnan(ice_mask)] = np.nan
        sit_pro[sit_pro == 1] = 0
        sit_pro[sit_pro == 2] = 0
        sit_pro[sit_pro == 3] = 1
        sit_pro[sit_pro == 4] = 0.5
        sit_pro[sit_pro == 255] = np.nan

    water = np.ma.array(np.zeros_like(nsidc_age), mask=nsidc_age != 0)
    #raise

    plt.close('all')
    # nsidc factor = 2
    r1,r2 = 150,540
    c1,c2 = 180,540
    # nsidc factor = 4
    r1,r2 = 300,1080
    c1,c2 = 360,1080

    if doplot: fig,axs = plt.subplots(2, 4)

    #if doplot: ax00 = axs[0,0].imshow(nsidc_age[r1:r2, c1:c2], cmap='jet', vmin=1, vmax=4)
    #if doplot: ax0wm = axs[0,0].imshow(water[r1:r2, c1:c2], cmap='gray')
    #if doplot: axs[0,0].set_title('SIA NSIDC', fontsize=8)
    #nsidc_age_list.append(nsidc_age[r1:r2, c1:c2])

    nsidc_myi = np.array(nsidc_age)
    nsidc_myi[nsidc_myi==1] = 0
    nsidc_myi[nsidc_myi > 0] = 1
    if doplot:
        ax00 = axs[0,0].imshow(nsidc_myi[r1:r2, c1:c2], cmap='jet', vmin=minc, vmax=maxc)
        ax0wm = axs[0,0].imshow(water[r1:r2, c1:c2], cmap='gray')
        axs[0,0].set_title('MYI NSIDC/NSIDC', fontsize=8)
    nsidc_age_list.append(nsidc_age[r1:r2, c1:c2])

    #agem_osi_pro_predom = np.argmax(agew_osi_pro, axis=0) + 1.
    #agem_osi_pro_predom[np.isnan(ice_mask)] = np.nan
    #if doplot: ax01 = axs[0,1].imshow(agem_osi_pro_predom[r1:r2, c1:c2], cmap='jet', vmin=1, vmax=4)
    #if doplot: ax2wm = axs[0,1].imshow(water[r1:r2, c1:c2], cmap='gray')
    #if doplot: axs[0,1].set_title('SIA/P NERSC', fontsize=8)
    #nersc_agep_list.append(agem_osi_pro_predom[r1:r2, c1:c2])

    myi_nsidc = agew_nsidc[1:].sum(axis=0)
    myi_nsidc[np.isnan(ice_mask)] = np.nan
    if doplot:
        ax01 = axs[0,1].imshow(myi_nsidc[r1:r2, c1:c2], cmap='jet', vmin=0, vmax=1)
        ax3wm = axs[0,1].imshow(water[r1:r2, c1:c2], cmap='gray')
        axs[0,1].set_title('MYI NSIDC/NERSC', fontsize=8)
    nsidc_agew_list.append(agew_nsidc[:, r1:r2, c1:c2])

    # show BREMEN MYI
    brem_file = glob.glob(os.path.join(brem_idir, nsidc_date.strftime('ECICE-ASCAT-AMSR2-%Y%m*.nc')))
    if len(brem_file) > 0:
        nmyi = Dataset(brem_file[0]).variables['nmyi'][:].reshape(896, 608)
        nmyi_pro = reproject_ice(brem_dom, nsidc_dom, nmyi) / 100.
        if doplot:
            ax02 = axs[0,2].imshow(nmyi_pro[r1:r2, c1:c2], cmap='jet', vmin=minc, vmax=maxc)
            ax2wm = axs[0,2].imshow(water[r1:r2, c1:c2], cmap='gray')
            axs[0,2].set_title('MYI UNI-BREMEN', fontsize=8)

    if doosisaf:
        if doplot:
            ax03 = axs[0,3].imshow(sit_pro[r1:r2, c1:c2], cmap='jet', vmin=minc, vmax=maxc)
            ax1wm = axs[0,3].imshow(water[r1:r2, c1:c2], cmap='gray')
            axs[0,3].set_title('MYI OSISAF', fontsize=8)
    #osisaf_sit_list.append(sit_pro[r1:r2, c1:c2])

    nersc_agew_list.append(agew_osi_pro[:, r1:r2, c1:c2])
    
    age1 = agew_osi_pro[0]
    age1[np.isnan(ice_mask)] = np.nan

    if len(agew_osi_pro) > 1:
        age2 = agew_osi_pro[1]
    else:
        age2 = np.zeros(ice_mask.shape)
    age2[np.isnan(ice_mask)] = np.nan
    if doplot:
        ax10 = axs[1,0].imshow(age2[r1:r2, c1:c2], cmap='jet', vmin=0, vmax=1)
        ax3wm = axs[1,0].imshow(water[r1:r2, c1:c2], cmap='gray')
        axs[1,0].set_title('2YO OSISAF/NERSC', fontsize=8)

    if len(agew_osi_pro) > 2:
        age3 = agew_osi_pro[2]
    else:
        age3 = np.zeros(ice_mask.shape)
    age3[np.isnan(ice_mask)] = np.nan
    if doplot:
        ax11 = axs[1,1].imshow(age3[r1:r2, c1:c2], cmap='jet', vmin=0, vmax=1)
        ax3wm = axs[1,1].imshow(water[r1:r2, c1:c2], cmap='gray')
        axs[1,1].set_title('3YO OSISAF/NERSC', fontsize=8)

    if len(agew_osi_pro) > 3:
        age4 = agew_osi_pro[3]
    else:
        age4 = np.zeros(ice_mask.shape)
    age4[np.isnan(ice_mask)] = np.nan
    if doplot:
        ax12 = axs[1,2].imshow(age4[r1:r2, c1:c2], cmap='jet', vmin=0, vmax=1)
        ax3wm = axs[1,2].imshow(water[r1:r2, c1:c2], cmap='gray')
        axs[1,2].set_title('4YO OSISAF/NERSC', fontsize=8)

    myi = agew_osi_pro[1:].sum(axis=0)
    myi[np.isnan(ice_mask)] = np.nan
    if doplot:
        ax13 = axs[1,3].imshow(myi[r1:r2, c1:c2], cmap='jet', vmin=0, vmax=1)
        ax3wm = axs[1,3].imshow(water[r1:r2, c1:c2], cmap='gray')
        axs[1,3].set_title('MYI OSISAF/NERSC', fontsize=8)

    if doplot: 
        [(ax.set_xticks([]),ax.set_yticks([])) for ax in axs.flatten()]

        #cbaxes = fig.add_axes([0.0, 0, 0.45, 0.03]) 
        #cb = fig.colorbar(ax00, cax=cbaxes, orientation='horizontal')
        #cb.ax.tick_params(labelsize=8)
        #cb.set_label('SIA, %s' % nsidc_date.strftime('%Y-%m-%d'), size=8)

        cbaxes2 = fig.add_axes([0.0, 0, 0.9, 0.03]) 
        cb2 = fig.colorbar(ax10, cax=cbaxes2, orientation='horizontal')
        cb2.ax.tick_params(labelsize=8)
        cb2.set_label('MYI fraction, %s' % nsidc_date.strftime('%Y-%m-%d'), size=8)

        plt.tight_layout()
        plt.savefig(odir + 'new_age_maps_%03d_conc.png' % k, dpi=200, bbox_inches='tight', pad_inches=0)
        plt.close()
    k += 1


labels = [f.split('.')[3][2:]+'/'+f.split('.')[4] for f in nsidc_files]

nsidc_age_list = np.array(nsidc_age_list)
nsidc_agew_list = np.array(nsidc_agew_list)
nersc_agew_list = np.array(nersc_agew_list)


AREA = 6.250 * 6.250 / 1000000.
XTICSTEP = 13
### NSIDC AGE
colors = ['#000000', '#00007f', '#00d4ff', '#f4db00', '#7f0000']

age_ts, age_name = (nsidc_age_list, 'NSIDC_NSIDC')
age_ts[nsidc_age_list==0]=0
cumsum = np.zeros(age_ts.shape[0])
x = range(age_ts.shape[0])
for age in range(0, 4):
    cumsum1 = cumsum + (age_ts == age).sum(axis=1).sum(axis=1)
    plt.fill_between(x, cumsum * AREA, cumsum1 * AREA, label='%sYI' % age, color=colors[age])
    cumsum = cumsum1

cumsum1 = cumsum + (age_ts > age).sum(axis=1).sum(axis=1)
plt.fill_between(x, cumsum * AREA, cumsum1 * AREA, label='%s+YI' % (age+1), color=colors[age+1])
cumsum = cumsum1

plt.ylim()
plt.legend(loc=3)
plt.title('Weekly %s SIA' % age_name)
plt.xticks(range(0,len(nsidc_files),XTICSTEP), labels[0:len(nsidc_files):XTICSTEP])
plt.xlabel('years/weeks')
plt.ylabel('area, $10^6$ km')

plt.savefig(odir + 'weekly_%s_SIA.png' % age_name, dpi=200, bbox_inches='tight', pad_inches=0)
plt.close()


AGES = 5
for agew_list, agew_name in [(nsidc_agew_list, 'NSIDC_NERSC'), (nersc_agew_list, 'OSI_NERSC')]:
    areas_list = np.zeros((AGES, len(agew_list)))
    age_areas_list = []
    for i, agew in enumerate(agew_list):
        print agew.shape
        agew[:, np.isnan(ice_mask[r1:r2, c1:c2])] = np.nan
        wat = 1. - agew.sum(axis=0)
        agew_wat = np.array([wat]+list(agew))
        areas = np.nansum(agew_wat, axis=(1,2))
        areas_list[:len(areas), i] = areas
    
    areas_list_sum = np.cumsum(areas_list, axis=0) * AREA
    for i in reversed(range(len(areas_list_sum))):
        plt.fill_between(x, areas_list_sum[i], 0, label='%sYI' % i, color=colors[i])
    plt.legend(loc=3)
    plt.title('Weekly %s SIA' % agew_name)
    plt.xticks(range(0,len(nsidc_files),XTICSTEP), labels[0:len(nsidc_files):XTICSTEP])
    plt.xlabel('years/weeks')
    plt.ylabel('area, $10^6$ km')
    plt.savefig(odir + 'weekly_%s_SIA.png' % agew_name, dpi=200, bbox_inches='tight', pad_inches=0)
    plt.close()
