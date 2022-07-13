import glob

import os
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from scipy.ndimage import zoom, distance_transform_edt, generic_filter, median_filter
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
        plt.imshow(replace_w_mean)
        plt.show()
        
        fig, ax = plt.subplots(1,2, figsize=(14,7))
        ax[0].imshow(dxo, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)
        ax[1].imshow(dyo, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)
        plt.show()
    
    return dxo, dyo
    
def smooth(dxo, dyo, sic_mask, showfig = False):
    '''Generic filter to remove outliers and blur nearest neighbors
    '''
    dxof = generic_filter(dxo, np.nanmedian, 3)
    dyof = generic_filter(dyo, np.nanmedian, 3)

    dxoff = fill_gaps_nn(dxof, 5)
    dyoff = fill_gaps_nn(dyof, 5)

    dxoff[np.isnan(sic_mask[1:-1:3, 1:-1:3])] = np.nan
    dyoff[np.isnan(sic_mask[1:-1:3, 1:-1:3])] = np.nan

    
    if showfig:
        fig, ax = plt.subplots(1,3,figsize=(21,7))
        ax[0].imshow(dyof, interpolation='nearest', clim=[-20, 20], cmap='gray')
        ax[0].imshow(sic_mask[1:-1:3, 1:-1:3], interpolation='nearest', clim=[0,1], cmap='bwr', alpha=0.2)

        ax[1].imshow(dxoff, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)
        ax[2].imshow(dyoff, interpolation='nearest', clim=[-20, 20], cmap=cm.balance)
        plt.show()

    return dxoff, dyoff
        
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
    

def get_files_for_mean(sid_files, index, year, month, sid_dir):
    '''return 3 sid files for average around the define date
    [index-1:index+1]
    if 1st or last day of the month,
        get the last day of previous month, or the 1st day of next month
    '''
    if index == 0:  # first day of 1 month
        pmo = get_previous_month(mo) 
        pyr = int(year) - 1
        if pyr < 1991:  # min year
            return sid_files[:3]
        pre_sid = [f for f in sorted(glob.glob(f'{sid_dir}/{pyr}/{pmo}/*/*nc'))][-1]      
        return [pre_sid] + sid_files[index:index+2]
    
    elif index == len(sid_files)-1:  # last day of 1 month
        nmo = get_next_month(mo)
        nyr = int(year) + 1
        if nyr > 2020:  # max year
            return sid_files[index-2:index+1]
        next_sid = [f for f in sorted(glob.glob(f'{sid_dir}/{nyr}/{nmo}/*/*nc'))][0]
        return sid_files[index-1:index+1] + [next_sid]

    else:  # all other cases
        return sid_files[index-1:index+2]


def get_previous_month(mo):
    '''return previous month in string format
    '''
    prev_month = int(mo) - 1
    if prev_month < 1:
        return '12'
    else:
        return f'{prev_month}'.zfill(2)

def get_next_month(mo):
    '''return next month in string format
    '''
    next_month = int(mo) + 1
    if next_month > 12:
        return '01'
    else:
        return f'{next_month}'.zfill(2)
    
# files location
sic_dir = '/Data/sim/data/OSISAF_ice_conc_CDR/'
sid_dir = '/Data/sim/data/OSISAF_ice_drift_CDR_v1pre_multi_oi/'
outdir = '/home/leoede/sea_ice_age/newosi_quality_test/'

# flags parameters
flag_values = np.array([ 0,  1,  2,  3,  4,  5, 10, 11, 12, 13, 20, 21, 22, 23, 30])
flag_conv = np.zeros(31)
flag_conv[flag_values] = np.arange(flag_values.size)


# loop on all years, all months
# years = [str(x) for x in range(1991,2021)]
years = [str(x) for x in range(2014,2018)]

months = [f'{x}'.zfill(2) for x in range(1,13)]

# years = ['1991','1992','1993']  # for testing
# months = ['03','04']

# create all subfolders to store results
isdir_or_create(years, outdir)
               
        
for year in years:
    for mo in months:
        sic_files = sorted(glob.glob(f'{sic_dir}/{year}/{mo}/*nc')) # [:3]
        sic_dates = [f.split('_')[-1].split('.')[0] for f in sic_files]
        sid_files = [f for f in sorted(glob.glob(f'{sid_dir}/{year}/{mo}/*/*nc')) if f.split('-')[-1].split('.')[0] in sic_dates]

        # define output dir
        odir = f'{outdir}{year}/'
        odir_png = f'{outdir}{year}/png/'

        # loop on all files 
        for idx, (sic_file, sid_file) in enumerate(zip(sic_files, sid_files)):
            # get files for mean on 3 days
            m_sid_files = get_files_for_mean(sid_files, idx, year, mo, sid_dir)
            # retrieve mean dx, dy
            ddx, ddy = retrieve_mean_dxdy(m_sid_files)
            # retrieve sic, dx, dy
            sic, dx, dy, dflag, sic_mask = retrieve_file_variables(sic_file, sid_file)
            # first correction
            dxo, dyo = replace_pixels(dx, dy, dflag, ddx, ddy, showfig=False)
            # second correction
            dxoff, dyoff = smooth(dxo, dyo, sic_mask)
            # save as .npz
            ofile = os.path.join(odir, os.path.basename(sid_file) + '.npz')
            print('Save:', os.path.basename(ofile))
            np.savez(ofile, u=dxoff, v=dyoff, c=sic)
            # plot png: u, v, sic
            plot_opngs(ofile, odir_png)
