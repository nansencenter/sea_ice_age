import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from ovl_plugins.lib.lagrangian import rungekutta4

#ftp://osisaf.met.no/archive_test/ice/drift_lr/osi405c_demo_archive/2013/01/ice_drift_nh_polstere-625_multi-oi_201301191200-201301211200.nc
#wget -w 1 -r -nc -nd -A 'ice_drift_nh_*.nc' -P /Data/sat/downloads/osi405c_demo_archive/ ftp://osisaf.met.no/archive_test/ice/drift_lr/osi405c_demo_archive/2016


#### IMPLEMENT THE NSIDC ICE AGE ALGORITHM
from iceagelib import *

idir_uv = '/files/nsidc0116_icemotion_vectors_v3/'
idir_ia = '/files/nsidc0611_seaice_age_v3/'
odir = '/files/sea_ice_age/nsidc_f02_xy/'

factor=2
res = 25000

#icemotion.grid.week.1979.07.n.v3.bin
ifiles = sorted(glob.glob(idir_uv + 'icemotion.grid.week*.bin'))[:321]

"""
propagate_from(get_i_of_file_nsidc(1978, 44, ifiles), ifiles, factor=factor, saveice=False, savexy=True, odir=odir)
propagate_from(get_i_of_file_nsidc(1979, 35, ifiles), ifiles, factor=factor, saveice=False, savexy=True, odir=odir)
propagate_from(get_i_of_file_nsidc(1980, 35, ifiles), ifiles, factor=factor, saveice=False, savexy=True, odir=odir)
propagate_from(get_i_of_file_nsidc(1981, 35, ifiles), ifiles, factor=factor, saveice=False, savexy=True, odir=odir)
propagate_from(get_i_of_file_nsidc(1982, 35, ifiles), ifiles, factor=factor, saveice=False, savexy=True, odir=odir)
propagate_from(get_i_of_file_nsidc(1983, 35, ifiles), ifiles, factor=factor, saveice=False, savexy=True, odir=odir)
propagate_from(get_i_of_file_nsidc(1984, 35, ifiles), ifiles, factor=factor, saveice=False, savexy=True, odir=odir)
"""
rfiles = sorted(glob.glob(odir + '*1984.52.n.v3.bin_icemap_xy.npz'), reverse=True)

ices = []
for rfile in rfiles:
    x1 = np.load(rfile)['x1']
    y1 = np.load(rfile)['y1']
    #plt.plot(x1.flat,y1.flat,'.');plt.show()

    c1 = np.floor(x1 * factor / res).astype('uint32')
    r1 = 722 - np.floor(y1 * factor / res).astype('uint32')

    gpi = (c1>0) * (c1<1000)
    #plt.plot(c1[gpi], r1[gpi],'.');plt.show()

    idx = c1[gpi] + r1[gpi] * r1.shape[1]
    bincount = np.bincount(idx)
    ice = np.zeros_like(c1)
    ice.flat[:bincount.size] = bincount
    ices.append(ice.astype(np.float32))

ices = np.array(ices)
ices_frac = np.divide(ices, ices.sum(axis=0)) * 100
ices_frac[np.isnan(ices_frac)] = 0

# load NSIDC age
idir_ia = '/files/nsidc0611_seaice_age_v3/'
nsidc_f = glob.glob(idir_ia + 'iceage.grid.week.1984.52.n.v3.bin')[0]
nsidc_age = np.fromfile(nsidc_f, np.uint8).reshape(361*2,361*2).astype(np.float32)
nsidc_age[nsidc_age==255] = np.nan
nsidc_age /= 5.

r1,r2 = 150,450
c1,c2 = 200,500

fig, axs = plt.subplots(2,4)
age = 2
img_axs = []
for ice_frac, ax in zip(ices_frac, axs.flat):
    ice_frac[np.isnan(nsidc_age)] = np.nan
    img_axs.append(ax.imshow(ice_frac[r1:r2,c1:c2], vmin=0, vmax=100, cmap='jet'))
    ax.set_title('%d-yo fraction' % age, fontsize=8)
    age += 1
    
img_axs.append(axs.flat[-1].imshow(nsidc_age[r1:r2,c1:c2], cmap='jet'))
axs.flat[-1].set_title('SIA, 1985.01.01', fontsize=8)

[(ax.set_xticks([]),ax.set_yticks([])) for ax in axs.flat]

cbaxes = fig.add_axes([0.0, 0.1, 0.45, 0.03]) 
cb = fig.colorbar(img_axs[0], cax=cbaxes, orientation='horizontal')
cb.ax.tick_params(labelsize=8)
cb.set_label('fraction, %', size=8)

cbaxes2 = fig.add_axes([0.5, 0.1, 0.45, 0.03]) 
cb2 = fig.colorbar(img_axs[-1], cax=cbaxes2, orientation='horizontal')
cb2.ax.tick_params(labelsize=8)
cb2.set_label('SIA, years', size=8)


plt.tight_layout()
plt.savefig('NSIDC_SIA_fractions.png', dpi=200, bbox_inches='tight', pad_inches=0)
plt.close()

