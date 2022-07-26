import glob

import os
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from scipy.ndimage import zoom, distance_transform_edt, generic_filter, median_filter, gaussian_filter
import cmocean.cm as cm

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
    dist, indi = distance_transform_edt(mask, return_distances=True, return_indices=True)
    gpi = dist <= distance
    r, c = indi[:, gpi]
    array = np.array(array)
    array[gpi] = array[r, c]
    return array

def gaussian_filter_nan(a, sigma, truncate):
    b = a.copy()
    b[np.isnan(a)]=0
    bb = gaussian_filter(b, sigma=sigma, truncate=truncate)
    w = np.ones_like(b)
    w[np.isnan(b)] = 0
    ww = gaussian_filter(w, sigma=sigma, truncate=truncate)
    return bb / ww


def retrieve_mean_dxdy(sid_files):
    '''
    returns ddx, dx mean of all sid files
            ddy, dy mean of all sid files
    '''
    ddx = []
    ddy = []

    for sid_file in sid_files:
        print(sid_file)
        dds = Dataset(sid_file)
        dx = dds['dX'][0].filled(np.nan)
        dy = dds['dY'][0].filled(np.nan)

        ddx.append(dx)
        ddy.append(dy)

    return np.nanmean(ddx, axis=0), np.nanmean(ddy, axis=0)


def retrieve_nanmean_dxdy(sid_files):
    '''2nd iteration of correction: we take the mean and we put nan to bad pixels values (detected with flag)
    '''
    ddx = []
    ddy = []
    ddflag = []
    
    ddx_nan = []
    ddy_nan = []

    # get data from files
    for sid_file in sid_files:
        print(sid_file)
        dds = Dataset(sid_file)
        ddx.append(dds['dX'][0].filled(np.nan))
        ddy.append(dds['dY'][0].filled(np.nan))
        ddflag.append(dds['status_flag'][0])
        
    # average with nan
    for dx,dy,dflag in zip(ddx, ddy, ddflag):
        bad_pix = (dflag==22) * (dx == 0) * (dy == 0)
        dx[bad_pix] = np.nan
        ddx_nan.append(dx)
        dy[bad_pix] = np.nan
        ddy_nan.append(dy)

    dx_mean = np.nanmean(ddx_nan, axis=0)
    dy_mean = np.nanmean(ddy_nan, axis=0)
    
    return dx_mean, dy_mean


def replace_pixels(dx, dy, dflag, dx_mean, dy_mean, showfig = False):
    '''Replace pixels with drift x and y == 0 and flag == 22
    Replace by nearest neighbors averaged with mean
    '''  
    replace_w_mean = (dflag==22) * (dx == 0) * (dy == 0)  # pixels to replace

    dxo = np.array(dx)
    dyo = np.array(dy)

    dx[np.isnan(dx)] = 0
    dy[np.isnan(dy)] = 0
    dx[replace_w_mean] = np.nan
    dy[replace_w_mean] = np.nan
    dx = fill_gaps_nn(dx, 100)
    dy = fill_gaps_nn(dy, 100)

    dx = (dx + dx_mean) / 2
    dy = (dy + dy_mean) / 2

    dxo[replace_w_mean] = dx[replace_w_mean]
    dyo[replace_w_mean] = dy[replace_w_mean]


    if showfig:
        fig, ax = plt.subplots(1,3, figsize=(21,7))
        ax[0].imshow(replace_w_mean)
        ax[1].imshow(dxo, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)
        ax[2].imshow(dyo, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)
        plt.savefig(f"/home/leoede/sea_ice_age/newosi_quality_test/test/replace_pixels.png")
        plt.close()
    
    return dxo, dyo, replace_w_mean
    
def smooth(dxo, dyo, sic_mask, savefig = False):
    '''Generic filter to remove outliers and blur nearest neighbors
    '''
    dxof = generic_filter(dxo, np.nanmedian, 3)
    dyof = generic_filter(dyo, np.nanmedian, 3)

    dxof = gaussian_filter_nan(dxof, 1, 1)
    dyof = gaussian_filter_nan(dyof, 1, 1)

    dxoff = fill_gaps_nn(dxof, 5)
    dyoff = fill_gaps_nn(dyof, 5)

    dxoff[np.isnan(sic_mask[1:-1:3, 1:-1:3])] = np.nan
    dyoff[np.isnan(sic_mask[1:-1:3, 1:-1:3])] = np.nan

    
    if savefig:
        fig, ax = plt.subplots(1,3,figsize=(21,7))
        ax[0].imshow(dyof, interpolation='nearest', clim=[-20, 20], cmap='gray')
        ax[0].imshow(sic_mask[1:-1:3, 1:-1:3], interpolation='nearest', clim=[0,1], cmap='bwr', alpha=0.2)

        ax[1].imshow(dxoff, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)
        ax[2].imshow(dyoff, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)
        plt.savefig(f"/home/leoede/sea_ice_age/newosi_quality_test/test/gen_smooth.png")

    return dxoff, dyoff

def advanced_smooth(dxo, dyo, sic_mask, showfig = False):
    '''Second generation of filter: gaussian filter
    '''
    dxof1 = gaussian_filter_nan(dxo, 1, 2)
    dyof1 = gaussian_filter_nan(dyo, 1, 2)

    dxoff1 = fill_gaps_nn(dxof1, 5)
    dyoff1 = fill_gaps_nn(dyof1, 5)

    dxoff1[np.isnan(sic_mask[1:-1:3, 1:-1:3])] = np.nan
    dyoff1[np.isnan(sic_mask[1:-1:3, 1:-1:3])] = np.nan
    
    if showfig:
        fig, ax = plt.subplots(1,3,figsize=(21,7))
        ax[0].imshow(dyof1, interpolation='nearest', clim=[-20, 20], cmap='gray')
        ax[0].imshow(sic_mask[1:-1:3, 1:-1:3], interpolation='nearest', clim=[0,1], cmap='bwr', alpha=0.2)

        ax[1].imshow(dxoff1, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)
        ax[2].imshow(dyoff1, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)
        plt.savefig(f"/home/leoede/sea_ice_age/newosi_quality_test/test/adv_smooth.png")
        plt.close()
    
    return dxoff1, dyoff1
    
def retrieve_file_variables(sic_file, sid_file):
    '''returns sic, dx, dy, dflag, sic_mask for one pair of file
    '''
    cds = Dataset(sic_file)
    dds = Dataset(sid_file)
    sic = cds['ice_conc'][0].filled(np.nan)
    dx = dds['dX'][0].filled(np.nan)
    dy = dds['dY'][0].filled(np.nan)
    dflag = dds['status_flag'][0]

    sic_mask = (sic > 0).astype(float)
    sic_mask[sic_mask == 0] = np.nan
    
    return sic, dx, dy, dflag, sic_mask
    
def plot_opngs(ifile, odir_png, umax=20):
    npz = np.load(ifile)
    ofile = f'{odir_png}/{os.path.basename(ifile)}.png'
    print('Plot', ifile)
    
    fig, ax = plt.subplots(1,3,figsize=(15,5))
    ax[0].imshow(npz['u'], clim=[-umax, umax], cmap=cm.balance, interpolation='nearest')
    ax[1].imshow(npz['v'], clim=[-umax, umax], cmap=cm.balance, interpolation='nearest')
    ax[2].imshow(npz['c'], interpolation='nearest')

    plt.tight_layout()
    plt.savefig(ofile, bbox_inches='tight', pad_inches=0, dpi=150)
    plt.close()    

def isdir_or_create(years, outdir):
    '''Check a subfolder exists for all year in years
    and check subsubfolder '/png' exists
    
    create the folders when non existant
    
    years, list of string
    '''
    for year in years:
        npz_dir = f'{outdir}{year}/'
        png_dir = f'{outdir}{year}/png/'
        if not os.path.isdir(npz_dir):
            os.makedirs(npz_dir)
        if not os.path.isdir(png_dir):
            os.makedirs(png_dir)

    return
    
def get_files_for_mean(sid_files, index, year, sid_dir):
    '''return 3 sid files for average around the define date (+/- 3 day)
    if 1st or last day of the month,
        get the last day of previous month, or the 1st day of next month
    '''
    if index == 0 and year == '1991':  # first data available: (+ 2 day)
        return sid_files[:7]   
    elif index == len(sid_files)-1 and year == '2020':  # last data available: (- 2 day)
        return sid_files[-7:]
    else:  # all other cases: 
#         return sid_files[index-1:index+2]  # (+/- 1 day)
        return sid_files[index-3:index+4]  # (+/- 1 day)

def compute_dist_weight(replace_w_mean, max_dist=10, showfig=False):
    '''
    '''
    dist = distance_transform_edt(replace_w_mean, return_distances=True) + distance_transform_edt(~replace_w_mean, return_distances=True)
    dist_weight = 1 - np.clip(dist, 0, max_dist)/max_dist
    
    if showfig:
        fig, ax = plt.subplots(1,1,figsize=(5,5))
        plt.imshow(dist_weight)
        plt.colorbar()
        plt.savefig(f"/home/leoede/sea_ice_age/newosi_quality_test/test/dist_weight.png")
        plt.close()
    
    return dist_weight

def apply_weight_correction(dist_weight, dxoff1, dyoff1, dxoff, dyoff, savefig=False):
    '''Apply distance weighted correction to dx, dy
    '''
    dxw = dist_weight * dxoff1 + (1 - dist_weight) * dxoff
    dyw = dist_weight * dyoff1 + (1 - dist_weight) * dyoff

    if savefig:
        fig, ax = plt.subplots(2,3,figsize=(21,15))
        ax[0,0].imshow(dxoff, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)
        ax[1,0].imshow(dyoff, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)

        ax[0,1].imshow(dxw, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)
        ax[1,1].imshow(dyw, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)

        ax[0,2].imshow(dxw-dxoff, interpolation='nearest', clim=[-5, 5], cmap=cm.balance)
        ax[1,2].imshow(dyw-dyoff, interpolation='nearest', clim=[-5, 5], cmap=cm.balance)
        plt.savefig(f"/home/leoede/sea_ice_age/newosi_quality_test/test/weight_correction.png")
        plt.close()
    
    return dxw, dyw

    

# files location
sic_dir = '/Data/sim/data/OSISAF_ice_conc_CDR/'
sid_dir = '/Data/sim/data/OSISAF_ice_drift_CDR_v1pre_multi_oi/'
outdir = '/Data/sim/data/OSISAF_ice_drift_CDR_v1pre_multi_oi/postproc_v3/'


# loop on all years, all months
years = [str(x) for x in range(1991,2021)]

# create all subfolders to store results
isdir_or_create(years, outdir)
               
sic_files = sorted(glob.glob(f'{sic_dir}/*/*/*nc'))# [:3]
# remove sic_files before 1991: because sid files not available
# remove after 2021 01 01 for same reason
beg_idx = np.where(np.array([sfile.find('1991') for sfile in sic_files]) >-1)[0][0]
end_idx = np.where(np.array([sfile.find('2021/') for sfile in sic_files]) >-1)[0][0]
sic_files = sic_files[beg_idx:end_idx]

sic_dates = [f.split('_')[-1].split('.')[0] for f in sic_files]
sid_files = [f for f in sorted(glob.glob(f'{sid_dir}/*/*/*/*nc')) if f.split('-')[-1].split('.')[0] in sic_dates]

# loop on all files 
for idx, (sic_file, sid_file) in enumerate(zip(sic_files, sid_files)):
    # get year and output directories
    year = sid_file[51:55]
    odir = f'{outdir}{year}/'
    odir_png = f'{outdir}{year}/png/'

    # if ofile exists, skip processing
    ofile = os.path.join(odir, os.path.basename(sid_file) + '.npz')
#     if os.path.isfile(ofile):
#         print(f'Skipping process for: {ofile}')
#         continue

    # get files for mean on 7 days
    m_sid_files = get_files_for_mean(sid_files, idx, year, sid_dir)
    # retrieve mean dx, dy
    ddx, ddy = retrieve_nanmean_dxdy(m_sid_files)
    # retrieve sic, dx, dy
    sic, dx, dy, dflag, sic_mask = retrieve_file_variables(sic_file, sid_file)
    # first correction
    dxo, dyo, replace_w_mean = replace_pixels(dx, dy, dflag, ddx, ddy)
    # second correction
    dxoff, dyoff = smooth(dxo, dyo, sic_mask)
    dxoff1, dyoff1 = advanced_smooth(dxo, dyo, sic_mask)
    # correction as a function of the distance to the bad pixels
    dist_weight = compute_dist_weight(replace_w_mean, max_dist=10)
    dxw, dyw = apply_weight_correction(dist_weight, dxoff1, dyoff1, dxoff, dyoff)
    
    # save as .npz
    print('Save:', os.path.basename(ofile))
#     np.savez(ofile, u=dxoff, v=dyoff, c=sic)  # first correction: not enough
    np.savez(ofile, u=dxw, v=dyw, c=sic)
    
    # plot png: u, v, sic
    plot_opngs(ofile, odir_png)
