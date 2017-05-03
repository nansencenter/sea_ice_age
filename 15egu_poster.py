import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from nansat import *

from iceagelib import *

idir = '/files/sea_ice_age/osi_newprop_f10_zoom1_conc/'

# NSIDC ICE AGE
nsidc_nsr = NSR('+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0 +no_defs')
nsidc_dom = Domain(nsidc_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 12500 12500')

# domain of OSISAF concentration
osi_conc_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_conc_dom = Domain(osi_conc_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 10000 10000')

# OSI SIE ICE AGE
osi_nsr='+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45'
osi_dom = Domain(osi_nsr, '-te -3781250 -5281250 3656250 5781250 -tr 6250 6250')

# DESTINATION DOMAIN    
dst_nsr = NSR('+proj=stere +datum=WGS84 +ellps=WGS84 +lat_0=90 +lon_0=0  +lat_ts=70 +no_defs')
#dst_dom = Domain(dst_nsr, '-te -4512500 -4512500 4512500 4512500 -tr 6250 6250')
dst_dom = Domain(dst_nsr, '-te -2325000 -2675000 2050000 2325000 -tr 6250 6250')

#"""
### PLOT ICE FRACTIONS
ifiles = [
    'icemap_2012-10-01_2012-10-02.npz',
    'icemap_2012-10-01_2013-03-15.npz',
    'icemap_2012-10-01_2013-09-17.npz',
    'icemap_2012-10-01_2014-03-17.npz',
    'icemap_2013-09-17_2013-09-18.npz',
    'icemap_2013-09-17_2014-03-17.npz',
]

plt.close('all')
ice_pros = []
for ifile in ifiles:
    print ifile
    ice = np.load(idir + ifile)['ice']
    ice_pro = reproject_ice(osi_dom, dst_dom, ice)
    #ice_pro[ice_pro > 0] = 1 - ice_pro[ice_pro > 0]
    #ice_pro[:] = 0
    nmap = Nansatmap(dst_dom, resolution='l')
    nmap.imshow(ice_pro, cmap='magma', vmin=0, vmax=1)
    nmap.save(ifile + '.png', dpi=250)
    ice_pros.append(ice_pro)
#"""

myi_2014_03_17 = ice_pros[3] + ice_pros[5]
fyi_2014_03_17 = myi_2014_03_17.copy()
fyi_2014_03_17[fyi_2014_03_17 > 0] = 1 - fyi_2014_03_17[fyi_2014_03_17 > 0]

nmap = Nansatmap(dst_dom, resolution='l')
nmap.imshow(myi_2014_03_17, cmap='magma', vmin=0, vmax=1)
nmap.save('icemap_myi_2014_03_17.png', dpi=250)

nmap = Nansatmap(dst_dom, resolution='l')
nmap.imshow(fyi_2014_03_17, cmap='magma', vmin=0, vmax=1)
nmap.save('icemap_fyi_2014_03_17.png', dpi=250)

'''
montage 2012-10-02_conc.png\
        2013-03-15_conc.png\
        2013-09-17_conc.png\
        2014-03-17_conc.png\
        -tile 4x1 -geometry +0+0 mont01.png 

montage icemap_2012-10-01_2012-10-02.npz.png\
         icemap_2012-10-01_2013-03-15.npz.png\
         icemap_2012-10-01_2013-09-17.npz.png\
         icemap_2012-10-01_2014-03-17.npz.png\
         -tile 4x1 -geometry +0+0 mont02.png 

montage icemap_2012-10-01_2012-10-02.npz_zero.png\
        icemap_2012-10-01_2013-03-15.npz_diff.png\
        icemap_2013-09-17_2013-09-18.npz.png\
        icemap_2013-09-17_2014-03-17.npz.png\
        -tile 4x1 -geometry +0+0 mont03.png 

montage icemap_2012-10-01_2012-10-02.npz_zero.png\
        icemap_2012-10-01_2012-10-02.npz_zero.png\
        icemap_2012-10-01_2012-10-02.npz_zero.png\
        icemap_fyi_2014_03_17.png\
        -tile 4x1 -geometry +0+0 mont04.png 

montage 2012-10-02_age.png\
        2013-03-15_age.png\
        2013-09-17_age.png\
        2014-03-17_age.png\
        -tile 4x1 -geometry +0+0 mont05.png 

'''
raise 
### PLOT ICE WEIGHTED AGE
idates = [
    dt.datetime(2012,10,02),
    dt.datetime(2013,3,15),
    dt.datetime(2013,9,17),
    dt.datetime(2014,3,17),
]

for idate in idates:
    agew_osi, agem_osi = get_mean_age(idir, idate)
    ice_conc = reproject_ice(osi_conc_dom, osi_dom, get_ice_conc(idate))
    ice_mask = ice_conc > 15
    agew_osi, agem_osi = add_fyi(agew_osi, agem_osi, ice_mask)
    agem_osi_pro = reproject_ice(osi_dom, dst_dom, agem_osi)
    ice_conc_pro = reproject_ice(osi_dom, dst_dom, ice_conc)

    nmap = Nansatmap(dst_dom, resolution='l')
    nmap.imshow(agem_osi_pro, cmap='jet', vmin=0, vmax=3)
    nmap.save(idate.strftime('%Y-%m-%d') + '_age.png', dpi=250)
    plt.close('all')

    nmap = Nansatmap(dst_dom, resolution='l')
    nmap.imshow(ice_conc_pro, cmap='magma', vmin=0, vmax=100)
    nmap.save(idate.strftime('%Y-%m-%d') + '_conc.png', dpi=250)
    plt.close('all')

plt.close('all')
fig = plt.figure(figsize=(4, 0.5))
ax1 = fig.add_axes([0.05, 0.5, 0.9, 0.4])
cmap = mpl.cm.magma
norm = mpl.colors.Normalize(vmin=0, vmax=100)
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.set_label('Sea ice fraction, %')
plt.savefig('fraction_colorbar.png', dpi=150)


plt.close()
fig = plt.figure(figsize=(4, 0.5))
ax1 = fig.add_axes([0.05, 0.5, 0.9, 0.4])
cmap = mpl.cm.jet
norm = mpl.colors.Normalize(vmin=0, vmax=3)
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.set_label('Sea ice age, years')
plt.savefig('age_colorbar.png', dpi=150)
plt.close()


idir_ia = '/files/nsidc0611_seaice_age_v3/'
nsidc_files = sorted(glob.glob(idir_ia + 'iceage.grid.week.2014.11.n.v3.bin'))
for nsidc_file in nsidc_files:
    nsidc_age = np.fromfile(nsidc_file, np.uint8).reshape(361*2,361*2).astype(np.float32)
    nsidc_age[nsidc_age==255] = np.nan
    nsidc_age /= 5.
    nsidc_age_pro = reproject_ice(nsidc_dom, dst_dom, nsidc_age)    

    nmap = Nansatmap(dst_dom, resolution='l')
    nmap.imshow(nsidc_age_pro, cmap='jet', vmin=0, vmax=3)
    nmap.save(os.path.basename(nsidc_file) + '.png', dpi=250)
    plt.close('all')
