import os
import glob
import datetime as dt

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

from nansat import Nansat
from ovl_plugins.lib.interpolation import fill_gaps_nn
from ovl_plugins.lib.lagrangian import rungekutta4
from scipy.ndimage.interpolation import zoom

#### IMPLEMENT THE NSIDC ICE AGE ALGORITHM

def zoom_nan(img, factor, order=1):
    imgz = np.array(img, np.float32)
    mask = np.isfinite(imgz).astype(np.float32)
    imgz = fill_gaps_nn(imgz, 3)
    imgz[np.isnan(imgz)] = 0
    
    imgz = zoom(imgz, factor, order=order)
    mask = zoom(mask, factor, order=order)
    imgz[mask < 0.5] = np.nan
    return imgz
    

def read_uv_nsidc(ifile, factor=1):
    d = np.fromfile(ifile, np.int16)
    u = d[0::3].reshape(361,361) / 1000. # m/s
    v = d[1::3].reshape(361,361) / 1000. # m/s
    c = d[2::3].reshape(361,361)
    c[c > 0] = 100
    u[c == 0] = np.nan
    v[c == 0] = np.nan

    if factor != 1:
        u = zoom_nan(u, factor)
        v = zoom_nan(v, factor)
        c = zoom_nan(c, factor)
        
    return u,v,c

def get_i_of_file_nsidc(year, week, ifiles):
    for i, ifile in enumerate(ifiles):
        if '%04d.%0d' % (year, week) in ifile:
            return i

def get_i_of_file_osi(year, month, day, ifiles):
    for i, ifile in enumerate(ifiles):
        if os.path.basename(ifile).split('_')[-1].startswith('%04d%02d%02d' % (year, month, day)):
            return i

def get_ice_conc(date, n=0):
    print date, n
    url_template = 'http://thredds.met.no/thredds/dodsC/myocean/siw-tac/siw-metno-glo-osisaf/conc/%Y/%m/ice_conc_nh_polstere-100_multi_%Y%m%d%H%S.nc'
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

def read_uv_osi_filled(ifile):
    u = np.load(ifile)['u']
    v = np.load(ifile)['v']
    c = np.load(ifile)['ice']

    return u,v,c
    
def reproject_ice(d0, d1, ice0):
    n = Nansat(domain=d0, array=ice0)
    n.reproject(d1, addmask=False, blockSize=10)
    return n[1]
    
    
def propagate_from(i0, ifiles, reader=read_uv_nsidc, res=25000, factor=2,
                    h=60*60*24*7, repro=None, odir='./'):
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

def get_nsidc_date(nsidc_file):
    y = int(os.path.basename(nsidc_file).split('.')[3])
    w = int(os.path.basename(nsidc_file).split('.')[4])
    return dt.datetime(y, 1, 1) + dt.timedelta(w*7)
    
def propagate_from_newprop(i0, ifiles, reader=read_uv_nsidc, res=25000, factor=1,
                    h=60*60*24*7, repro=None, odir='./'):

    u, v, f = read_uv_nsidc(ifiles[i0], factor=factor)

    y0 = ifiles[i0].split('.')[3]
    w0 = ifiles[i0].split('.')[4]
    ice0files = sorted(glob.glob(odir + '*%s.%s.n.v3.bin_icemap.npy.npz' % (y0, w0)))

    ice0 = np.zeros(f.shape)  # water is zero age old
    ice0[f > 0] = 1 # on 1978, week 44 all ice is in second year

    #import ipdb; ipdb.set_trace()

    # reduce 100% concentration by previous concentrations
    for ice0file in ice0files:
        y0prev = ice0file.split('.')[3]
        if y0prev < y0:
            print 'Correct ICE0 for', y0prev, ice0file
            ice0prev = np.load(ice0file)['ice']
            ice0 -= ice0prev

    ice0[f <= 0] = 0
    #ice0[ice0 < 0] = 0

    # save initiation day YYYY.WW-YYYY.WW
    ofile = '%s/%s_%s_icemap' % (odir,
                        os.path.basename(ifiles[i0]),
                        os.path.basename(ifiles[i0]))
    np.savez_compressed(ofile+'.npy', ice=ice0)
    
    cols0, rows0 = np.meshgrid(range(u.shape[1]), range(u.shape[0]))
    k = 0
    for i in range(i0, len(ifiles)-1):
        print os.path.basename(ifiles[i0]), os.path.basename(ifiles[i])
        u, v, f = reader(ifiles[i], factor=factor)

        dc = u * h / res
        dr = - v * h / res

        cols1 = cols0+dc
        rows1 = rows0+dr

        dca = 1 - (cols1 - np.floor(cols1))
        dcb = cols1 - np.floor(cols1)
        dra = 1 - (rows1 - np.floor(rows1))
        drb = rows1 - np.floor(rows1)

        cols1a = (np.floor(cols1)).astype(np.uint32)
        cols1b = (np.floor(cols1) + 1).astype(np.uint32)
        rows1a = (np.floor(rows1)).astype(np.uint32)
        rows1b = (np.floor(rows1) + 1).astype(np.uint32)

        saa = dca * dra
        sab = dca * drb
        sba = dcb * dra
        sbb = dcb * drb

        ice1aa = saa * ice0
        ice1ab = sab * ice0
        ice1ba = sba * ice0
        ice1bb = sbb * ice0

        gpi = np.isfinite(cols1) * np.isfinite(rows1)
        ice1 = np.zeros_like(ice0)

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
        
        next_year = ifiles[i+1].split('.')[3]
        next_week = ifiles[i+1].split('.')[4]
        prev_files = sorted(glob.glob(odir + '*%s.%s.n.v3.bin_icemap.npy.npz' % (next_year, next_week)))
        sum_prev_ice1 = np.zeros_like(ice1)

        print 'Sum previous: ', 
        for prev_file in prev_files:
            prev_date = get_nsidc_date(prev_file)
            if prev_date < get_nsidc_date(ifiles[i0]):
                print prev_date, 
                sum_prev_ice1 += np.load(prev_file)['ice']
        print
        max_ice_increment = 1-sum_prev_ice1

        #if k > 15:
        #    import ipdb; ipdb.set_trace()

        ice1[ice1 > max_ice_increment] = max_ice_increment[ice1 > max_ice_increment]
        
        ice0 = ice1

        ofile = '%s/%s_%s_icemap' % (odir,
                            os.path.basename(ifiles[i0]),
                            os.path.basename(ifiles[i+1]))
        np.savez_compressed(ofile+'.npy', ice=ice0)
        
        k += 1

def collect_age(rfiles):
    age0 = np.load(rfiles[0])['ice']
    ice_age = np.ones(age0.shape)
    theage = 2
    for rfile in rfiles:
        age = np.load(rfile)['ice']
        ice_age[age > 0] = theage
        theage += 1
    
    return ice_age

def vis_ice_npz(pref, vmin=0, vmax=1):
    ifiles = sorted(glob.glob(pref+'*.npz'))
    k = 0
    for ifile in ifiles:
        ice = np.load(ifile)['ice']
        plt.imsave('%s_frame_%05d.png' % (pref, k), ice, cmap='jet', vmin=vmin, vmax=vmax)
        k += 1

def vis_drift_npz(idir, odir):
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
