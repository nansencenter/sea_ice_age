import os
import glob
import datetime as dt
from dateutil.parser import parse

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from netCDF4 import Dataset
import pygrib

from nansat import Nansat, Nansatmap
from ovl_plugins.lib.interpolation import fill_gaps_nn
from ovl_plugins.lib.lagrangian import rungekutta4
from scipy.ndimage.interpolation import zoom
from scipy.ndimage.filters import gaussian_filter, maximum_filter


#### IMPLEMENT THE NSIDC ICE AGE ALGORITHM
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

def zoom_nan(img, factor, order=1, sigma=0):
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
    imgz = fill_gaps_nn(imgz, 10)
    imgz[np.isnan(imgz)] = 0
    
    imgz = zoom(imgz, factor, order=order)
    mask = zoom(mask, factor, order=order)
    imgz[mask < 0.5] = np.nan
    
    if sigma > 0:
        imgz = nangaussian_filter(imgz, sigma)
    
    return imgz
    
def get_nsidc_uv(ifile, src_res=25000, h=7*24*60*60, factor=1, order=1, sigma=0, **kwargs):
    ''' Load U,V and ice mask from NSIDC file
    Return U,V in row/col pixels
    '''
    d = np.fromfile(ifile, np.int16)
    u = d[0::3].reshape(361,361) / 1000. # m/s
    v = d[1::3].reshape(361,361) / 1000. # m/s
    c = d[2::3].reshape(361,361)
    c[c > 0] = 100
    u[c == 0] = np.nan
    v[c == 0] = np.nan
    
    if factor != 1:
        u = zoom_nan(u, factor, order=order, sigma=sigma)
        v = zoom_nan(v, factor, order=order, sigma=sigma)
        c = zoom_nan(c, factor, order=order, sigma=sigma)

    u =  u * factor * h / src_res # pix
    v = -v * factor * h / src_res # pix

    return u, v, c

def get_nsidc_date(nsidc_file):
    ''' Get date of NSIDC file '''
    y = int(os.path.basename(nsidc_file).split('.')[3])
    w = int(os.path.basename(nsidc_file).split('.')[4])
    return dt.datetime(y, 1, 1) + dt.timedelta(w*7)

def get_nsidc_i_of_file(year, week, ifiles):
    ''' Find index of NSIDC file based on input date'''
    for i, ifile in enumerate(ifiles):
        if '%04d.%02d' % (year, week) in ifile:
            return i

def get_osi_i_of_file(year, month, day, ifiles):
    ''' Find index of OSISAF file based on input date'''
    for i, ifile in enumerate(ifiles):
        if os.path.basename(ifile).split('_')[-1].startswith('%04d%02d%02d' % (year, month, day)):
            return i

def get_osi_sic(idir, idate, n=0):
    ''' Get sea ice concentration from downloaded GRIB file '''
    sic_files = glob.glob(idir + 'ice_conc*%s????.grb' % idate.strftime('%Y%m%d'))

    if len(sic_files) == 0:
        print 'NO ', idate, 
        idate += dt.timedelta(1)
        print 'TRY ', idate
        return get_osi_sic(idir, idate, n=n+1)

    ifile = sic_files[0]
    grbs = pygrib.open(ifile)
    for grb in grbs:
        if 'Ice cover' in grb['name']:
            break
    c = np.flipud(grb.data()[0])
    if grbs.messagenumber == 1:
        # OLD
        pole_mask = c.mask
        pole_mask = maximum_filter(pole_mask, 3)
        c_f = fill_gaps_nn(c.data, 10, pole_mask)
        c_f[np.isnan(c_f)] = 100
        c_f[c_f < 0] = np.nan
    else:
        # NEW
        c_f = c.data
        c_f[c.mask] = np.nan
        
    return c_f

def get_osi_uvf(ifile):
    ''' Load U,V status flag from OSISAF LR drift file '''
    n0 = Nansat(ifile, mapperName='generic')
    status_flag = n0['status_flag']
    u = n0['dX']*1000/60/60/48
    if n0.has_band('dY_v1p4'):
        dY = n0['dY_v1p4']
    else:
        dY = n0['dY']
    v = dY*1000/60/60/48
    return u, v, status_flag

def fill_osi_uv(u, v, c, sid_dom, sic_dom, nn_dst=5, sigma=2, **kwargs):
    ''' Resample and fill gaps in UV data using C'''
    uv_pro = []
    for uv in [u,v]:
        # upscale U,V to the grid of C
        uvp = reproject_ice(sid_dom, sic_dom, uv)
        # fill LAND with 0 speed
        uvp[np.isnan(c)] = 0
        # extrapolate U,V (from ice and land onto empty pixels)
        uvp = fill_gaps_nn(uvp, nn_dst)
        # replace WATER pixels with NAN
        uvp[c == 0] = np.nan
        # replace LAND pixels with NAN
        uvp[np.isnan(c)] = np.nan
        # smooth upscaled U,V
        uvp = nangaussian_filter(uvp, sigma)
        uv_pro.append(uvp)

    return uv_pro[0], uv_pro[1], c

def get_osi_uvc_filled(sid_file, src_res=10000, h=24*60*60, factor=1, **kwargs):
    ''' Load presaved SID and SIC '''
    data = np.load(sid_file)
    u, v, c = data['u'], data['v'], data['c']
    u =  u * factor * h / src_res # pix
    v = -v * factor * h / src_res # pix

    return u, v, c

def get_osi_date(osi_file):
    ''' Get date of OSISAF file '''
    return parse(os.path.basename(osi_file).split('_')[-1][:8])
    
def reproject_ice(d0, d1, ice0, eResampleAlg=0):
    ''' Convert ice product from one projection to another '''
    n = Nansat(domain=d0, array=ice0)
    n.reproject(d1, addmask=False, eResampleAlg=eResampleAlg, blockSize=10)
    return n[1]

def propagate_fowler(i_start, i_end, ifiles, reader, get_date, odir,
                     savexy=False, **kwargs):
    ''' Apply NSIDC algorithm for ice age 
    Input:
        i0, index of file to start from
        ifiles: list of files with U,V,C
        reader: function to read U,V,C
        get_date: function to read date
        src_res: spatial resolution of source dataset (m)
        factor: zoom factor
        h: time step (sec)
        odir: output directory
    Output:
        None. Saves files with sea ice age.
    '''
    rnd = lambda x: np.round(x).astype(int)
    # read initial conditions
    u, v, c = reader(ifiles[i_start], **kwargs)
    d0 = get_date(ifiles[i_start])

    # save initiation day YYYY.WW-YYYY.WW
    ice0 = np.zeros(c.shape) + np.nan
    ice0[c > 0] = 1
    ofile = '%s/icemap_%s_%s.npz' % (odir,
                        d0.strftime('%Y-%m-%d'),
                        d0.strftime('%Y-%m-%d'))
    np.savez_compressed(ofile, ice=ice0)

    # define initial coordinates
    cols, rows = np.meshgrid(np.arange(u.shape[1], dtype=float),
                             np.arange(u.shape[0], dtype=float))
    gpi = np.isfinite(cols) # all pixels are good
    # loop through input files
    for ifile in ifiles[i_start+1:i_end+1]:
        print 'PROP: ', os.path.basename(ifile)
        cols_gpi = cols[gpi] + u[rnd(rows[gpi]), rnd(cols[gpi])]
        rows_gpi = rows[gpi] + v[rnd(rows[gpi]), rnd(cols[gpi])]
        cols[gpi] = cols_gpi
        rows[gpi] = rows_gpi
        cols[~gpi] = np.nan
        rows[~gpi] = np.nan
        
        gpi = (np.isfinite(cols * rows) *
               (cols >= 0) *
               (rows >= 0) *
               (cols < u.shape[1]) *
               (rows < u.shape[0]))

        ice = np.zeros(u.shape) + np.nan
        u, v, c = reader(ifile, **kwargs)
        d = get_date(ifile)
        ice[c > 0] = 0        
        ice[rnd(rows[gpi]), rnd(cols[gpi])] = 1
        ofile = '%s/icemap_%s_%s.npz' % (odir,
                            d0.strftime('%Y-%m-%d'),
                            d.strftime('%Y-%m-%d'))
        np.savez_compressed(ofile, ice=ice)
        if savexy:
            np.savez_compressed(ofile+'_xy.npz', rows=rows, cols=cols)

def get_icemap_dates(icemap_file):
    ''' Get dates of from ice map filename''' 
    nameparts = os.path.splitext(os.path.basename(icemap_file))[0].split('_')
    d0 = parse(nameparts[-2])
    d1 = parse(nameparts[-1])
    return d0, d1

def propagate_nersc(i_start, i_end, ifiles, reader, get_date, odir, 
                    conc=False, min_flux=0, **kwargs):
    ''' Apply NERSC algorithm for ice age 
    Input:
        i_start, index of file to start from
        i_end, index of the last file
        ifiles: list of files with U,V,C
        reader: function to read U,V,C from input file
        get_date: function to read date from input file
        src_res: spatial resolution (m)
        factor: zoom factor
        h: time step (sec)
        odir: output directory
    Output:
        None. Files with sea ice age.
    '''
    # get initial ice mask and drift
    u0, v0, c0 = reader(ifiles[i_start], **kwargs)
    #import ipdb; ipdb.set_trace()
    d0 = get_date(ifiles[i_start])
    
    # find files with ice age produced from previous years
    ice0files = sorted(glob.glob(odir + '*%s.npz' % d0.strftime('%Y-%m-%d')))

    # initialize ice age fraction map and set all concentrations to 100% 
    if conc:
        ice0 = c0 / 100.
    else:
        ice0 = np.zeros(c0.shape) # water is zero years old
        ice0[c0 > 0] = 1 # on 15 September all ice is 1 year old
    
    # reduce initial concentration by concentrations of older ice
    for ice0file in ice0files:
        ice0d0, ice0d1 = get_icemap_dates(ice0file)
        if ice0d0.year < d0.year:
            print 'Correct ICE0 for', ice0d0.year
            ice0prev = np.load(ice0file)['ice']
            ice0 -= ice0prev

    # fix negative concentration
    #ice0[f0 <= 0] = 0
    ice0[ice0 < 0] = 0

    # save initiation day YYYY.WW-YYYY.WW
    ofile = '%s/icemap_%s_%s.npz' % (odir,
                        d0.strftime('%Y-%m-%d'),
                        d0.strftime('%Y-%m-%d'))
    np.savez_compressed(ofile, ice=ice0)
    
    # initial coordinates of each pixel
    cols0, rows0 = np.meshgrid(range(u0.shape[1]), range(u0.shape[0]))
    k = 0
    # loop over files with ice drift
    for i in range(i_start, i_end):
        d1 = get_date(ifiles[i+1])
        ofile = '%s/icemap_%s_%s.npz' % (odir,
                            d0.strftime('%Y-%m-%d'),
                            d1.strftime('%Y-%m-%d'))
        # skip processed
        if os.path.exists(ofile):
            ice0 = np.load(ofile)['ice']
            continue
            
        print os.path.basename(ifiles[i_start]), os.path.basename(ifiles[i])
        # read U,V,C,T
        dc, dr, c = reader(ifiles[i], **kwargs)
        d = get_date(ifiles[i])

        # increment of coordinates
        #dc = u * h / res
        #dr = - v * h / res

        # allow only substantial enough increment
        #dc[np.abs(dc) < min_incr] = min_incr * np.sign(dc[np.abs(dc) < min_incr])
        #dr[np.abs(dr) < min_incr] = min_incr * np.sign(dr[np.abs(dr) < min_incr])
        
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
        
        # allow only substantial enough flux
        #for ice1flux in [ice1aa, ice1ab, ice1ba, ice1bb]:
        #    ice1flux[ice1flux < min_flux] = 0

        # find valid donor pixels
        gpi = np.isfinite(cols1) * np.isfinite(rows1)
        ice1 = np.zeros_like(ice0)

        # for four directions
        # w: ice flux from donor pixels
        # flat index of recipient pixels
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
            ice_conc = c / 100.
        else:
            ice_conc = np.ones_like(sum_prev_ice1)

        # maximum possible increment in sea ice fractional concentration
        max_ice_increment = ice_conc - sum_prev_ice1
        # limit increment by the maximum possible
        ice1[ice1 > max_ice_increment] = max_ice_increment[ice1 > max_ice_increment]
        
        ice0 = ice1
        np.savez_compressed(ofile, ice=ice0)

def collect_age(rfiles):
    ''' Compute MAX ice age from several ice age fractions ''' 
    age0 = np.load(rfiles[0])['ice']
    ice_age = age0 + 1 # FYI  = 1; MYI = 2
    theage = 2
    for rfile in rfiles:
        age = np.load(rfile)['ice']
        ice_age[age > 0] = theage
        theage += 1
    
    return ice_age

def vis_ice_npz(pref, vmin=0, vmax=1, nfiles=None):
    ''' Visualize NPZ file with single ice age fraction as a PNG '''
    ifiles = sorted(glob.glob(pref+'*.npz'))[:nfiles]
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

def get_mean_age(idir, thedate, ice_mask):
    ''' Compute weighted average of fractional ice age '''
    rfiles = sorted(glob.glob(idir + '*%s.npz' % thedate.strftime('%Y-%m-%d')), reverse=True)
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
    
    myi_weight = ice_age_weights.sum(axis=0)
    fyi_weight = 1 - myi_weight
    fyi_weight[ice_mask == 0] = 0
    
    ice_age_weights = np.array([fyi_weight] + list(ice_age_weights))
    ice_age_weighted = np.multiply(range(1,1+len(ice_age_weights)), ice_age_weights.T).T
    ice_age_mean = ice_age_weighted.sum(axis=0) / ice_age_weights.sum(axis=0)
    
    ice_age_mean[ice_mask == 0] = 0
    ice_age_mean[np.isnan(ice_mask)] = np.nan
    
    myi = np.sum(ice_age_weights[1:], axis=0)
    
    return ice_age_weights, ice_age_mean, myi

def save_max_age(idir, sia_factor, vmin=0, vmax=8):
    ''' Use Fowler method to compute SIA (max age in bin) '''
    odir = idir + 'sia/'
    if not os.path.exists(odir):
        os.makedirs(odir)

    ifiles = sorted(glob.glob(idir + '*.npz'))

    idates = map(get_icemap_dates, ifiles)
    dates0, dates1 = zip(*idates)
    dst_dates = sorted(set(dates1))

    k = 0
    for dst_date in dst_dates:
        print dst_date
        msk = '%s*%s' % (idir, dst_date.strftime('%Y-%m-%d.npz'))
        rfiles = sorted(glob.glob(msk), reverse=True)
        sia = collect_age(rfiles)

        if sia_factor == 2:
            ## REDUCE RESOLUTION OF NERSC PRODUCT 2 times
            sia = np.nanmax(np.stack([sia[0::2, 0::2],
                                      sia[1::2, 1::2]]), axis=0)
        elif sia_factor == 4:
            ## REDUCE RESOLUTION OF NERSC PRODUCT 4 times
            sia = np.nanmax(np.stack([sia[0::4, 0::4],
                                      sia[1::4, 1::4],
                                      sia[2::4, 2::4],
                                      sia[3::4, 3::4]]), axis=0)
        elif sia_factor == 5:
            ## REDUCE RESOLUTION OF NERSC PRODUCT 5 times
            sia = np.nanmax(np.stack([sia[0::5, 0::5],
                                      sia[1::5, 1::5],
                                      sia[2::5, 2::5],
                                      sia[3::5, 3::5],
                                      sia[4::5, 4::5]]), axis=0)

        ofile = '%s/%s_sia.npz' % (odir, dst_date.strftime('%Y-%m-%d'))
        np.savez_compressed(ofile, sia=sia)
        plt.imsave('%s/sia_%05d.png' % (odir, k), sia, cmap='jet', vmin=vmin, vmax=vmax)
        k += 1

def save_mean_age(sid_files, icemap_dir, reader, get_date, vmin=0, vmax=5, **kwargs):
    ''' Use NERSC method to compute SIA and ice fractions '''
    odir = icemap_dir + 'sia/'
    if not os.path.exists(odir):
        os.makedirs(odir)

    k = 0
    for sid_file in sid_files:
        u,v,c = reader(sid_file, **kwargs)
        d = get_date(sid_file)
        print d
        sif, sia, myi = get_mean_age(icemap_dir, d, c)
        ofile = '%s/%s_sia.npz' % (odir, d.strftime('%Y-%m-%d'))
        np.savez_compressed(ofile, sia=sia, sif=sif, myi=myi)
        plt.imsave('%s/sia_%05d.png' % (odir, k), sia, cmap='jet', vmin=vmin, vmax=vmax)
        k += 1
                    
def make_map(ifile, prod, src_dom, dst_dom, array=None,
             vmin=0, vmax=5, dpi=250, cmap='jet',
             title=None, text=None, textx=100000, texty=5000000, fontsize=18):
    ''' Make map with a product '''
    if array is None:
        sia = np.load(ifile)[prod]
    else:
        sia = array
    sia_pro = reproject_ice(src_dom, dst_dom, sia)
    nmap = Nansatmap(dst_dom, resolution='l')
    nmap.imshow(sia_pro, vmin=vmin, vmax=vmax, cmap=cmap)
    if title is not None:
        plt.title(title, fontsize=fontsize)
    if text is not None:
        plt.text(textx, texty, text, fontsize=fontsize, va='top', bbox=dict(facecolor='white', alpha=0.9))
    nmap.save('%s_%s.png' % (ifile, prod), dpi=dpi)

def save_legend(cmap, bounds, label, filename, format='%1i'):
    # colorbar for SIA
    cmap = plt.get_cmap(cmap)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    fig = plt.figure(figsize=(8, 1))
    ax2 = fig.add_axes([0.05, 0.5, 0.9, 0.3])
    cb = mpl.colorbar.ColorbarBase(ax2,
        cmap=cmap, norm=norm, spacing='proportional', ticks=bounds,
        boundaries=bounds, format=format, orientation='horizontal')
    cb.set_label(label, size=12)
    plt.savefig(filename, dpi=150, bbox_inches='tight', pad_inches=0)
    plt.close()

def get_nsidc_raw_sia(ifile):
    nsidc_age = np.fromfile(ifile, np.uint8).reshape(361*2,361*2).astype(np.float32)
    nsidc_age[nsidc_age==255] = np.nan
    nsidc_age /= 5.

    return nsidc_age
