import datetime as dt
import glob
import os

from dateutil.parser import parse
import matplotlib.pyplot as plt
from nansat import Nansat
import numpy as np
import pygrib
from scipy import ndimage
from scipy.ndimage.filters import generic_filter
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.filters import median_filter
from scipy.ndimage.interpolation import zoom


def get_osi_date(osi_file):
    ''' Get date of OSISAF file '''
    return parse(os.path.basename(osi_file).split('_')[-1][:8])

def nangaussian_filter(img, sigma):
    ''' Gaussian filter for images with NaN '''
    img0 = img.copy()
    img0[np.isnan(img)] = 0
    img0 = gaussian_filter(img0, sigma=sigma)

    imgw = 0 * img.copy() + 1
    imgw[np.isnan(imgw)] = 0
    imgw = gaussian_filter(imgw, sigma=sigma)
    
    img0 /= imgw
    img0[np.isnan(img)] = np.nan
    return img0

def fill_gaps_nn(array, distance=2, mask=None):
    """ Fill gaps in input array
    # https://stackoverflow.com/questions/5551286/filling-gaps-in-a-numpy-array/9262129#9262129
    Parameters
    ----------
    array : 2D numpy.array
        Raster with data
    distance : int
        Minimum size of gap to fill
    mask : 2D numpy.array
        Raster with True where to fill gaps
    Returns
    -------
    array : 2D numpy.array
        Raster with data with gaps filled
    """
    if mask is None:
        mask = np.isnan(array)
    dist, indi = ndimage.distance_transform_edt(mask, return_distances=True, return_indices=True)
    gpi = dist <= distance
    r, c = indi[:, gpi]
    array = np.array(array)
    array[gpi] = array[r, c]
    return array

def zoom_nan(img, factor, order=1, sigma=0):
    ''' Increase resolution of image with gaps  filled by nan
    Input:
        img: matrix to resize
        factor: resize factor
        order: interpolation order
        sigma: for gaussian filtering
    Output:
        img: resized image
    '''
    imgz = np.array(img, np.float32)
    mask = np.isfinite(imgz).astype(np.float32)
    imgz = fill_gaps_nn(imgz, 10)
    imgz[np.isnan(imgz)] = 0
    
    imgz = zoom(imgz, factor, order=order)
    mask = zoom(mask, factor, order=order)
    imgz[mask < 0.5] = np.nan
    
    if sigma > 0:
        imgz = nangaussian_filter(imgz, sigma)
    
    return imgz

def reproject_ice(d0, d1, ice0, resample_alg=0):
    ''' Convert ice product from one projection to another '''
    n = Nansat.from_domain(d0, array=ice0)
    n.reproject(d1, addmask=False, resample_alg=resample_alg, block_size=10)
    return n[1]

def fill_med_osi_uv_2(u, v, c, sid_dom_x, sic_dom, uvf0=None, zf=2, footprint=np.ones((2,7,7)), min_conc=0):
    ''' Median filter, resample and fill gaps in UV data using C and previous UV'''
    uv_pro = []
    uvf1 = []
    for i, uv in enumerate([u,v]):
        # increase resolution
        uvz = zoom_nan(uv, zf)
        # fill gaps with nearest neigbour
        uvf = fill_gaps_nn(uvz, 100)
        uvf1.append(uvf)
        # median filter (also previous UV, if exists)
        if uvf0 is None:
            uvm = median_filter(uvf, footprint=footprint[0], mode='nearest')
        else:
            uvm = median_filter([uvf0[i], uvf], footprint=footprint, mode='nearest')[1]
        
        # upscale U,V to the grid of C
        uvp = reproject_ice(sid_dom_x, sic_dom, uvm, 2)
        #import ipdb; ipdb.set_trace()
        # replace WATER pixels with NAN
        uvp[c < min_conc] = np.nan
        # replace LAND pixels with NAN
        uvp[np.isnan(c)] = np.nan
        uv_pro.append(uvp)

    return uv_pro[0], uv_pro[1], c, uvf1

def get_osi_sic(idir, idate, force_grb=False, n=0, improve=True):
    ''' Get sea ice concentration from downloaded GRIB file '''
    mask = f'{idir}/ice_conc*{idate.strftime("%Y%m%d")}????.grb'
    sic_files = glob.glob(mask)

    if len(sic_files) == 0:
        print('NO ', idate)
        idate += dt.timedelta(1)
        print('TRY ', idate)
        return get_osi_sic(idir, idate, force_grb=force_grb, n=n+1)

    ifile = sic_files[0]
    # try to load C from presaved (corrected) NPZ file
    ifile_npz = ifile + '.npz'
    if not force_grb and os.path.exists(ifile_npz):
        return np.load(ifile_npz)['c']

    grbs = pygrib.open(ifile)
    for grb in grbs:
        if 'Ice cover' in grb['name']:
            break
    c = np.flipud(grb.data()[0])
    c_f = c.data / 100.    
    if improve:
        if grbs.messagenumber == 1:
            # OLD
            pole_mask = c.mask
            pole_mask = maximum_filter(pole_mask, 3)
            c_f = fill_gaps_nn(c_f, 100, pole_mask)
            c_f[np.isnan(c_f)] = 1.
            c_f[c_f < 0] = np.nan
        else:
            # NEW
            c_f[c.mask] = np.nan
        
    return c_f

def fill_gap_sic(sic_dir, dateA, dateB, dateC):
    mask = f'{sic_dir}/ice_conc*{dateC}????.grb'
    fileC = glob.glob(mask)[0]
    fileC_npz = fileC + '.npz'
    if not os.path.exists(fileC_npz):
        cA = get_osi_sic(sic_dir, parse(str(dateA)), force_grb=True)
        cB = get_osi_sic(sic_dir, parse(str(dateB)), force_grb=True)
        
        cC = ((cA + cB) / 2)
        np.savez_compressed(fileC_npz, c=cC)

def get_osi_uvf(ifile):
    ''' Load U,V status flag from OSISAF LR drift file '''
    n0 = Nansat(ifile, mapperName='generic')
    status_flag = n0['status_flag']
    dX = n0['dX']
    if n0.has_band('dY_v1p4'):
        dY = n0['dY_v1p4']
    else:
        dY = n0['dY']
    u = dX*1000./48./60./60.
    v = dY*1000./48./60./60.
    return u, v, status_flag

def plot_pngs(ifile, odir, umax=0.15):
    npz = np.load(ifile)
    ofile = f'{odir}/{os.path.basename(ifile)}.png'
    print('Plot', ofile)
    landmask = np.isnan(npz['c']).astype(float)
    landmask[landmask == 0] = np.nan

    fig, ax = plt.subplots(1,3,figsize=(15,5))
    ax[0].imshow(npz['u'], clim=[-umax, umax], cmap='bwr', interpolation='nearest')
    ax[1].imshow(npz['v'], clim=[-umax, umax], cmap='bwr', interpolation='nearest')
    ax[2].imshow(npz['c'], interpolation='nearest')
    for a in ax:
        a.set_xlim([100, 600])
        a.set_ylim([800, 350])
        a.imshow(landmask, cmap='gray', alpha=0.5, interpolation='nearest')
    plt.tight_layout()
    plt.savefig(ofile, bbox_inches='tight', pad_inches=0, dpi=150)
    plt.close()