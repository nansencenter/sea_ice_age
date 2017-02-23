import os
import glob
import datetime as dt
from dateutil.parser import parse

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

from nansat import Nansat
from ovl_plugins.lib.interpolation import fill_gaps_nn
from ovl_plugins.lib.lagrangian import rungekutta4
from scipy.ndimage.interpolation import zoom

#### IMPLEMENT THE NSIDC ICE AGE ALGORITHM

def zoom_nan(img, factor, order=1):
    ''' Increase resolution of image with gaps  filled by nan
    Input:
        img: matrix to resize
        factor: resize factor
        order: interpolation order
    Output:
        img: resized image
    '''
    imgz = np.array(img, np.float32)
    mask = np.isfinite(imgz).astype(np.float32)
    imgz = fill_gaps_nn(imgz, 3)
    imgz[np.isnan(imgz)] = 0
    
    imgz = zoom(imgz, factor, order=order)
    mask = zoom(mask, factor, order=order)
    imgz[mask < 0.5] = np.nan
    return imgz
    

def read_uv_nsidc(ifile, factor=1, order=1):
    ''' Load U,V and ice mask from NSIDC file '''
    d = np.fromfile(ifile, np.int16)
    u = d[0::3].reshape(361,361) / 1000. # m/s
    v = d[1::3].reshape(361,361) / 1000. # m/s
    c = d[2::3].reshape(361,361)
    c[c > 0] = 100
    u[c == 0] = np.nan
    v[c == 0] = np.nan
    
    if factor != 1:
        u = zoom_nan(u, factor, order=order)
        v = zoom_nan(v, factor, order=order)
        c = zoom_nan(c, factor, order=order)

    return u, v, c

def get_nsidc_date(nsidc_file):
    ''' Get date of NSIDC file '''
    y = int(os.path.basename(nsidc_file).split('.')[3])
    w = int(os.path.basename(nsidc_file).split('.')[4])
    return dt.datetime(y, 1, 1) + dt.timedelta(w*7)

def get_i_of_file_nsidc(year, week, ifiles):
    ''' Find index of NSIDC file based on input date'''
    for i, ifile in enumerate(ifiles):
        if '%04d.%0d' % (year, week) in ifile:
            return i

def get_i_of_file_osi(year, month, day, ifiles):
    ''' Find index of OSISAF file based on input date'''
    for i, ifile in enumerate(ifiles):
        if os.path.basename(ifile).split('_')[-1].startswith('%04d%02d%02d' % (year, month, day)):
            return i

def get_ice_conc(date, n=0):
    ''' Get ice concentration from OSISAF for given date '''
    print date, n
    url_template = 'http://thredds.met.no/thredds/dodsC/myocean/siw-tac/siw-metno-glo-osisaf/conc/%Y/%m/ice_conc_nh_polstere-100_multi_%Y%m%d1200.nc'
    url = date.strftime(url_template)
    try:
        ds = Dataset(url)
    except IOError:
        date = date + dt.timedelta(1)
        return get_ice_conc(date, n=n+1)
    ice_conc = np.array(ds.variables['ice_conc'][:][0])
    status_flag = np.array(ds.variables['status_flag'][:][0])
    ds.close()
    ice_conc[status_flag == 100] = -100 # land
    ice_conc[status_flag == 101] = -101 # missing
    ice_conc[status_flag == 102] = -102 # unclassified (edge)
    return ice_conc
    

def read_uvc_osi(ifile, concDomain):
    ''' Load U,V and ice mask from OSISAF file '''
    print ifile
    n0 = Nansat(ifile)
    status_flag = n0['status_flag']
    u = n0['dX']*1000/60/60/48
    v = n0['dY']*1000/60/60/48

    datestr = os.path.basename(ifile).split('_')[-1].split('-')[0]
    date = dt.datetime(2010,1,1).strptime(datestr, '%Y%m%d%H%S')
    ice_conc = get_ice_conc(date)
    concN = Nansat(domain=concDomain, array=ice_conc)
    concN.reproject(n0, addmask=False)
    c = concN[1]
    
    return u,v,c

def read_uv_osi_filled(ifile, factor=1, order=1):
    ''' Load U,V and ice mask from OSISAF file with gaps filled '''
    u = np.load(ifile)['u']
    v = np.load(ifile)['v']
    c = np.load(ifile)['ice']
    u[c == 0] = np.nan
    v[c == 0] = np.nan

    if factor != 1:
        u = zoom_nan(u, factor, order=order)
        v = zoom_nan(v, factor, order=order)
        c = zoom_nan(c, factor, order=order)

    return u,v,c

def get_osi_date(osi_file):
    ''' Get date of OSISAF file '''
    return parse(os.path.basename(osi_file).split('_')[-1][:8])
    
def reproject_ice(d0, d1, ice0):
    ''' Convert ice product from one projection to another '''
    n = Nansat(domain=d0, array=ice0)
    n.reproject(d1, addmask=False, blockSize=10)
    return n[1]

def propagate_from(i0, ifiles, reader=read_uv_nsidc, res=25000, factor=2,
                    h=60*60*24*7, repro=None, odir='./'):
    ''' Apply NSIDC algorithm for ice age 
    Input:
        i0, index of file to start from
        ifiles: list of files with U,V,C
        reader: function to read U,V,C
        res: spatial resolution (m)
        factor: zoom factor
        h: time step (sec)
        repro: tuple with Domains for ice reprojection
        odir: output directory
    Output:
        None. Files with sea ice age.
    '''
        
    u0, v0, c0 = reader(ifiles[i0])

    x, y = np.meshgrid(range(u0.shape[1]), range(u0.shape[0], 0, -1))
    x *= res
    y *= res

    pt_cols, pt_rows = np.meshgrid(range(int(x.shape[1]*factor)),
                                   range(int(x.shape[0]*factor), 0, -1))
    x0 = pt_cols * res / factor
    y0 = pt_rows * res / factor
    # filter out water pixels
    #x0 = x0[c0 > 15]
    #y0 = y0[c0 > 15]

    for i in range(i0, len(ifiles)):
        print os.path.basename(ifiles[i0]), os.path.basename(ifiles[i])
        u1, v1, c1 = reader(ifiles[i])
        x1, y1 = rungekutta4(x, y, u0, v0, u1, v1, x0, y0, h)
        u0, v0, x0, y0 = u1, v1, x1, y1

        # U/V res
        c1 = (x1 / (res / factor)).astype(np.int16)
        r1 = (x.shape[0] * factor - y1 / (res / factor)).astype(np.int16)
        gpi = (np.isfinite(c1 * r1) *
               (c1 >= 0) *
               (r1 >= 0) *
               (c1 < x0.shape[1]) *
               (r1 < x0.shape[0]))
        ice1 = np.zeros(x0.shape)
        
        ice1[r1[gpi], c1[gpi]] = 1
        
        if repro is not None:
            ice1 = reproject_ice(repro[0], repro[1], ice1)

        ofile = '%s/%s_%s_icemap' % (odir,
                                    os.path.basename(ifiles[i0]),
                                    os.path.basename(ifiles[i]))
        np.savez_compressed(ofile+'.npy', ice=ice1.astype(np.bool))

def get_icemap_dates(icemap_file):
    ''' Get dates of from ice map filename''' 
    nameparts = os.path.splitext(os.path.basename(icemap_file))[0].split('_')
    d0 = parse(nameparts[-2])
    d1 = parse(nameparts[-1])
    return d0, d1

def propagate_from_newprop(i0, ifiles, reader=read_uv_nsidc, get_date=get_nsidc_date,
                           res=25000, factor=1, order=1,
                           h=60*60*24*7, repro=None, odir='./', conc=False):
    ''' Apply NERSC algorithm for ice age 
    Input:
        i0, index of file to start from
        ifiles: list of files with U,V,C
        reader: function to read U,V,C from input file
        get_date: function to read date from input file
        res: spatial resolution (m)
        factor: zoom factor
        h: time step (sec)
        repro: tuple with Domains for ice reprojection
        odir: output directory
    Output:
        None. Files with sea ice age.
    '''
    # get initial ice mask and drift
    u0, v0, f0 = reader(ifiles[i0], factor=factor)
    d0 = get_date(ifiles[i0])
    
    # find files with ice age produced from previous years
    ice0files = sorted(glob.glob(odir + '*%s.npz' % d0.strftime('%Y-%m-%d')))

    # initialize ice age fraction map and set all concentrations to 100% 
    ice0 = np.zeros(f0.shape)  # water is zero age old
    ice0[f0 > 0] = 1 # on 1 october all ice is 1 year old

    # reduce 100% concentration by concentrations of older ice
    for ice0file in ice0files:
        ice0d0, ice0d1 = get_icemap_dates(ice0file)
        if ice0d0.year < d0.year:
            print 'Correct ICE0 for', ice0d0.year
            ice0prev = np.load(ice0file)['ice']
            ice0 -= ice0prev

    # fix negative concentration
    ice0[f0 <= 0] = 0
    #ice0[ice0 < 0] = 0

    # save initiation day YYYY.WW-YYYY.WW
    ofile = '%s/icemap_%s_%s.npz' % (odir,
                        d0.strftime('%Y-%m-%d'),
                        d0.strftime('%Y-%m-%d'))
    np.savez_compressed(ofile, ice=ice0)
    
    # initial coordinates of each pixel
    cols0, rows0 = np.meshgrid(range(u0.shape[1]), range(u0.shape[0]))
    k = 0
    # loop over files with ice drift
    for i in range(i0, len(ifiles)-1):
        print os.path.basename(ifiles[i0]), os.path.basename(ifiles[i])
        # read U,V,C,T
        u, v, f = reader(ifiles[i], factor=factor)
        d = get_date(ifiles[i])

        # increment of coordinates
        dc = u * h / res
        dr = - v * h / res

        # coordinates at step 2
        cols1 = cols0+dc
        rows1 = rows0+dr

        # relative displacement of pixel
        dca = 1 - (cols1 - np.floor(cols1))
        dcb = cols1 - np.floor(cols1)
        dra = 1 - (rows1 - np.floor(rows1))
        drb = rows1 - np.floor(rows1)

        # corrdinates of recipient pixels
        cols1a = (np.floor(cols1)).astype(np.uint32)
        cols1b = (np.floor(cols1) + 1).astype(np.uint32)
        rows1a = (np.floor(rows1)).astype(np.uint32)
        rows1b = (np.floor(rows1) + 1).astype(np.uint32)

        # weights of ice flux into four directions
        # aa: 0
        # ab: 180
        # ba: 90
        # bb: -90
        saa = dca * dra
        sab = dca * drb
        sba = dcb * dra
        sbb = dcb * drb

        # flux into four directions
        ice1aa = saa * ice0
        ice1ab = sab * ice0
        ice1ba = sba * ice0
        ice1bb = sbb * ice0

        # find valid donor pixels
        gpi = np.isfinite(cols1) * np.isfinite(rows1)
        ice1 = np.zeros_like(ice0)

        # for four directions
        # w: ice flux from donor pixels
        # indeces of recipient pixels
        # bincount: sum of fluxes from donor pixels for given direction
        # ice1: sum of fluxes from donor pixels for all directions
        
        w =  ice1aa[rows0[gpi], cols0[gpi]]
        idx = cols1a[gpi] + rows1a[gpi] * rows1a.shape[1]
        bincount = np.bincount(idx, w)
        ice1.flat[:bincount.size] += bincount

        w =  ice1ab[rows0[gpi], cols0[gpi]]
        idx = cols1a[gpi] + rows1b[gpi] * rows1a.shape[1]
        bincount = np.bincount(idx, w)
        ice1.flat[:bincount.size] += bincount

        w =  ice1ba[rows0[gpi], cols0[gpi]]
        idx = cols1b[gpi] + rows1a[gpi] * rows1a.shape[1]
        bincount = np.bincount(idx, w)
        ice1.flat[:bincount.size] += bincount

        w =  ice1bb[rows0[gpi], cols0[gpi]]
        idx = cols1b[gpi] + rows1b[gpi] * rows1a.shape[1]
        bincount = np.bincount(idx, w)
        ice1.flat[:bincount.size] += bincount
        
        d1 = get_date(ifiles[i+1])
        
        # ice1: sum of fluxes for a given date
        # it should be corrected by older ice

        print 'Sum older: ', 
        # sum of fractions of older ice 
        sum_prev_ice1 = np.zeros_like(ice1)
        prev_files = sorted(glob.glob(odir + '*%s.npz' % d1.strftime('%Y-%m-%d')))
        # loop through files of older ice for the current date
        for prev_file in prev_files:
            prev_date, _ = get_icemap_dates(prev_file)
            if prev_date < d0:
                print prev_date, 
                sum_prev_ice1 += np.load(prev_file)['ice']
        print 'OK'
        
        # use SIC
        if conc:
            ice_conc = f / 100.
        else:
            ice_conc = np.ones_like(sum_prev_ice1)

        # maximum possible increment in sea ice fractional concentration
        max_ice_increment = ice_conc - sum_prev_ice1
        # limit increment by the maximum possible
        ice1[ice1 > max_ice_increment] = max_ice_increment[ice1 > max_ice_increment]

        ice0 = ice1

        ofile = '%s/icemap_%s_%s.npz' % (odir,
                            d0.strftime('%Y-%m-%d'),
                            d1.strftime('%Y-%m-%d'))
        np.savez_compressed(ofile, ice=ice0)

        k += 1

def collect_age(rfiles):
    ''' Compute MAX ice age from several ice age fractions ''' 
    age0 = np.load(rfiles[0])['ice']
    ice_age = np.ones(age0.shape)
    theage = 2
    for rfile in rfiles:
        age = np.load(rfile)['ice']
        ice_age[age > 0] = theage
        theage += 1
    
    return ice_age

def vis_ice_npz(pref, vmin=0, vmax=1):
    ''' Visualize NPZ file with single ice age fraction as a PNG '''
    ifiles = sorted(glob.glob(pref+'*.npz'))
    k = 0
    for ifile in ifiles:
        ice = np.load(ifile)['ice']
        plt.imsave('%s_frame_%05d.png' % (pref, k), ice, cmap='jet', vmin=vmin, vmax=vmax)
        k += 1

def vis_drift_npz(idir, odir):
    ''' Visualize drift '''
    ifiles = sorted(glob.glob('%s/*.npz' % idir))
    k = 0
    for ifile in ifiles:
        print ifile
        u =  np.load(ifile)['u'][::-1, :]
        v = -np.load(ifile)['v'][::-1, :]
        plt.quiver(u,v)
        plt.savefig('%s/uv_frame_%05d.png' % (odir, k), dpi=75, bbox_inches='tight', pad_inches=0)
        plt.close()
        k += 1

def get_mean_age(idir, thedate):
    ''' Compute weighted average of fractional ice age '''
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
    ''' Compute weight of FYI based on ice mask and MYI fractions '''
    myi_weight = ice_age_weights.sum(axis=0)
    fyi_weight = 1 - myi_weight
    fyi_weight[ice_mask == 0] = 0
    
    ice_age_weights = np.array([fyi_weight] + list(ice_age_weights))
    ice_age_weighted = np.multiply(range(1,1+len(ice_age_weights)), ice_age_weights.T).T
    ice_age_mean = ice_age_weighted.sum(axis=0) / ice_age_weights.sum(axis=0)
    
    ice_age_mean[ice_mask == 0] = 0
    ice_age_mean[np.isnan(ice_mask)] = np.nan
    
    return ice_age_weights, ice_age_mean
