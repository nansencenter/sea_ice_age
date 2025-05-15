import os
import glob

import numpy as np
from scipy.interpolate import RectBivariateSpline
from skimage.util import view_as_windows
from scipy.ndimage import uniform_filter

from remeshing import get_area
from utils import get_mesh_files, fill_gaps

def std_filter_using_windows(image, size=25):
    """
    Compute the standard deviation of an image using sliding windows.
    """
    pad_width = size // 2  # for 5x5 window
    image_pad = np.pad(image, pad_width, mode='constant', constant_values=np.nan)

    # Create a sliding window view of the image
    windows = view_as_windows(image_pad, (size, size), step=1)
    
    # Compute the standard deviation for each window
    std_image = np.nanstd(windows, axis=(2, 3))
    
    return std_image

def compute_c_std(ifile):
    # load cached gridded STD
    ofile = ifile.replace('ice_drift_nh', 'ice_conc_std_nh')
    if os.path.exists(ofile):
        c_std_grd = np.load(ofile, allow_pickle=True)['c_std'].item()
        return c_std_grd

    # compute new gridded STD
    c_grd = np.load(ifile)['c']
    c_grd[c_grd < 5] = np.nan
    c_std_grd = {}
    std_sizes = list(range(3,31,3))
    for size in std_sizes:
        c_std_grd[size] = std_filter_using_windows(c_grd[1:-1:3, 1:-1:3], size=size).astype(np.float16)
        c_std_grd[size][np.isnan(c_std_grd[size])] = 0
    
    # save to cache
    np.savez(ofile, c_std=c_std_grd)
    return c_std_grd 

class ComputeSidSicUncertainty:
    def __init__(self, ifiles, idates, start_date, sia_dir, mesh_dir, mesh_init_file, xc, yc):
        self.ifiles = ifiles
        self.idates = idates
        self.start_date = start_date
        self.sia_dir = sia_dir
        self.mesh_dir = mesh_dir
        self.mesh_init_file = mesh_init_file
        self.xc = xc
        self.yc = yc

    def compute_sid_sic_uncertainty(self, i):
        ifile = self.ifiles[i]
        idate = self.idates[i]

        ofile = f'{self.sia_dir}/unc/{self.start_date.year}/unc_sidsic_{self.start_date.strftime("%Y%m%d")}_{idate.strftime("%Y%m%d")}.npz'
        if os.path.exists(ofile):
            #unc_sic_sid = np.load(ofile)['unc_sic_sid']
            return #unc_sic_sid

        mesh_file, _ = get_mesh_files(idate, self.mesh_dir, self.mesh_init_file)
        with np.load(mesh_file) as d:
            x0 = d['x']
            y0 = d['y']
            t0 = d['t']
            x0el = x0[t0].mean(axis=1)
            y0el = y0[t0].mean(axis=1)

        unc_sid_file = f'{self.sia_dir}/unc/{self.start_date.year}/unc_sid_{self.start_date.strftime("%Y%m%d")}_{idate.strftime("%Y%m%d")}.npz'
        unc_sid_sum = np.load(unc_sid_file)['unc_sid']

        # compute uncertainty of SIC field as STD in the area defined by the uncert of drift field
        el_size = get_area(x0, y0, t0)**0.5
        el_size_mean = el_size.mean()
        sid_unc_rel_to_el_size = unc_sid_sum / el_size_mean

        min_rel_size_factor = 3
        min_rel_sizes = np.arange(3,31,3)-2

        c_std_grd = compute_c_std(ifile)
        unc_sic_sid = np.zeros_like(unc_sid_sum)
        for min_rel_size in min_rel_sizes:
            uncert_elems = np.nonzero(
                (unc_sid_sum > 0) *
                (sid_unc_rel_to_el_size >= min_rel_size) * 
                (sid_unc_rel_to_el_size < min_rel_size+min_rel_size_factor)
            )[0]

            x_unc_sid = x0el[uncert_elems]
            y_unc_sid = y0el[uncert_elems]

            grd_shape = self.yc[1:-1:3].size, self.xc[1:-1:3].size
            rbs = RectBivariateSpline(self.xc[1:-1:3], self.yc[1:-1:3], c_std_grd[min_rel_size+2][:grd_shape[0], :grd_shape[1]][::-1], kx=1, ky=1)
            c_std_msh = rbs(y_unc_sid, x_unc_sid, grid=False)  
            unc_sic_sid[uncert_elems] = c_std_msh

        # save to cache
        np.savez(ofile, unc_sic_sid=unc_sic_sid)
        #return unc_sic_sid

class ComputeSicUncertainty:
    def __init__(self, ifiles, idates, n_steps, mesh_init_file, mesh_dir, unc_dir, xc, yc):
        self.ifiles = ifiles
        self.idates = idates
        self.n_steps = n_steps
        self.mesh_init_file = mesh_init_file
        self.mesh_dir = mesh_dir
        self.unc_dir = unc_dir
        self.xc = xc
        self.yc = yc

    def load_sic_data(self, ifile, mesh_src_file):
        with np.load(ifile) as d:
            sic_unc = d['sic_unc']
        with np.load(mesh_src_file) as d:
            x = d['x']
            y = d['y']
            t = d['t']
            src2dst = d['src2dst']
            weights = d['weights']
        return sic_unc, x, y, t, src2dst, weights

    def compute_uncertainty(self, start_idx):
        stop_idx = min(start_idx + self.n_steps, len(self.idates))
        unc_sic_sum = None
        start_date = self.idates[start_idx]
        odir = f'{self.unc_dir}/{start_date.year}'
        print(f'Computing uncertainty for {start_date.strftime("%Y%m%d")}, {self.idates[stop_idx-1]}')
        os.makedirs(odir, exist_ok=True)
        for i in range(start_idx, stop_idx):
            ifile = self.ifiles[i]
            idate = self.idates[i]
            ofile = f'{odir}/unc_sic_{start_date.strftime("%Y%m%d")}_{idate.strftime("%Y%m%d")}.npz'
            if os.path.exists(ofile):
                with np.load(ofile) as d:
                    unc_sic_sum = d['unc_sic']
                    sic_min = d['sic_min']
                continue

            mesh_file, mesh_dst_file = get_mesh_files(idate, self.mesh_dir, self.mesh_init_file)
            unc_sic_grd, x, y, t, src2dst, weights = self.load_sic_data(ifile, mesh_file)
            unc_sic_fil = fill_gaps(unc_sic_grd, np.isnan(unc_sic_grd), distance=144)
            
            # interpolate SIC uncertainty to mesh
            unc_sic = RectBivariateSpline(self.xc, self.yc, unc_sic_fil[::-1], kx=1, ky=1)(y[t].mean(axis=1), x[t].mean(axis=1), grid=False)
            sic_file = mesh_file.replace('mesh', 'sic')
            c = np.load(sic_file)['c']

            if unc_sic_sum is None:
                # first time step
                unc_sic_sum = np.array(unc_sic)
                sic_min = np.array(c)
            else:
                # advect unc_sic_sum (from previous time step)
                unc_sic_sum_pro = np.zeros(src2dst[:,1].max()+1)
                np.add.at(unc_sic_sum_pro, src2dst[:,1], unc_sic_sum[src2dst[:,0]] * weights)
                
                # advect sic_min from previous time step
                sic_min_pro = np.zeros(src2dst[:,1].max()+1)
                np.add.at(sic_min_pro, src2dst[:,1], sic_min[src2dst[:,0]] * weights)

                # keep uncertainty of the minimal value of SIC
                min_sic_ids = c < sic_min_pro
                unc_sic_sum = unc_sic_sum_pro
                unc_sic_sum[min_sic_ids] = unc_sic[min_sic_ids]
                
                # keep minimal ice concentration
                sic_min = np.min([sic_min_pro, c], axis=0)
            np.savez(ofile, unc_sic=unc_sic_sum.astype(np.float16), sic_min=sic_min.astype(np.float16))

        #print(start_date)
        #fig, axs = plt.subplots(1, 1, figsize=(6, 6))
        #trp0 = axs.tripcolor(x, y, t, unc_sic_sum, cmap='jet')
        #plt.colorbar(trp0, ax=axs)
        #plt.show()            


class ComputeSidUncertainty:
    def __init__(self, ifiles, idates, mesh_dir, unc_dir, mesh_init_file, n_steps, size, xc, yc):
        self.ifiles = ifiles
        self.idates = idates
        self.mesh_dir = mesh_dir
        self.unc_dir = unc_dir
        self.mesh_init_file = mesh_init_file
        self.n_steps = n_steps
        self.size = size
        self.xc = xc
        self.yc = yc

    def load_sid_data(self, ifile, mesh_src_file):
        with np.load(ifile) as d:
            u = d['u']
            v = -d['v']
            sid_unc = d['sid_unc']
        u[np.isnan(u)] = 0
        v[np.isnan(v)] = 0
        with np.load(mesh_src_file) as d:
            x0 = d['x']
            y0 = d['y']
            t0 = d['t']
        return u, v, sid_unc, x0, y0, t0
    
    def compute_uncertainty(self, start_idx):
        stop_idx = min(start_idx + self.n_steps, len(self.idates))
        unc_sid_sum = None
        x = None
        start_date = self.idates[start_idx]
        odir = f'{self.unc_dir}/{start_date.year}'
        os.makedirs(odir, exist_ok=True)
        print('Processing', start_date, self.idates[stop_idx-1])

        for i in range(start_idx, stop_idx):
            ifile = self.ifiles[i]
            idate = self.idates[i]
            ofile = f'{odir}/unc_sid_{start_date.strftime("%Y%m%d")}_{idate.strftime("%Y%m%d")}.npz'
            if os.path.exists(ofile):
                unc_sid_sum = np.load(ofile)['unc_sid']
                continue

            mesh_file, _ = get_mesh_files(idate, self.mesh_dir, self.mesh_init_file)
            u_grd, _, unc_sid_grd, x, y, t = self.load_sid_data(ifile, mesh_file)
            unc_sid_fil = fill_gaps(unc_sid_grd, np.isnan(unc_sid_grd), distance=144)

            # compute uncertainty of drift field after smoothing
            unc_sid_smt = uniform_filter(unc_sid_fil**2, size=self.size) / self.size
            unc_sid_smt[u_grd == 0] = 0

            # interpolate uncertainty to mesh
            unc_sid = RectBivariateSpline(self.xc[1:-1:3], self.yc[1:-1:3], unc_sid_smt[::-1], kx=1, ky=1)(y[t].mean(axis=1), x[t].mean(axis=1), grid=False)

            if unc_sid_sum is None:
                # first time step
                unc_sid_sum = np.array(unc_sid)
            else:
                # advect unc_sid_sum (from previous time step)
                with np.load(mesh_file) as d:
                    src2dst = d['src2dst']
                    weights = d['weights']
                unc_sid_sum_pro = np.zeros(src2dst[:,1].max()+1)
                np.add.at(unc_sid_sum_pro, src2dst[:,1], unc_sid_sum[src2dst[:,0]] * weights)
                # compute uncertainty of drift field after advection
                unc_sid_sum = np.hypot(unc_sid, unc_sid_sum_pro)
            np.savez(ofile, unc_sid=unc_sid_sum.astype(np.float16))

class ComputeTotUncertainty:
    def __init__(self, mesh_dir, unc_dir, sid_dir, age_dir, mesh_init_file):
        self.mesh_dir = mesh_dir
        self.unc_dir = unc_dir
        self.sid_dir = sid_dir
        self.age_dir = age_dir
        self.mesh_init_file = mesh_init_file

        self.xc = np.load(self.mesh_init_file)['xc']
        self.yc = np.load(self.mesh_init_file)['yc']

    def load_age(self, idate):
        age_file = f'{self.age_dir}/{idate.year}/age_{idate.year}{idate.month:02d}{idate.day:02d}.npz'
        with np.load(age_file) as ds:
            a = ds['a']
            c = ds['c']
            f = ds['f']
            x = ds['x']
            y = ds['y']
            t = ds['t'].astype(int)
        return a, c, f, x, y, t

    def load_unc_obs(self, idate):
        ifile = sorted(glob.glob(f'{self.sid_dir}/{idate.year}/ice_drift*{idate.strftime("%Y%m%d")}*.npz'))[0]
        unc_sic_grd = np.load(ifile)['sic_unc']
        unc_sic_fil = fill_gaps(unc_sic_grd, np.isnan(unc_sic_grd), 100)
        unc_sic_fil[np.isnan(unc_sic_fil)] = 0
        a, c, f, x, y, t = self.load_age(idate)
        unc_obs = RectBivariateSpline(self.xc, self.yc, unc_sic_fil[::-1], kx=1, ky=1)(y[t].mean(axis=1), x[t].mean(axis=1), grid=False)
        return unc_obs

    def get_unc_sic_files(self, idate):
        unc_sic_files = []
        for i in range(7):
            src_year = idate.year - i
            unc_sic_file = f'{self.unc_dir}/{src_year}/unc_sic_{src_year}0905_{idate.year}{idate.month:02d}{idate.day:02d}.npz'
            if os.path.exists(unc_sic_file):
                unc_sic_files.append(unc_sic_file)
        return unc_sic_files

    def load_unc_myi(self, idate):
        unc_sic_files = self.get_unc_sic_files(idate)    
        unc_myi = []
        for unc_sic_file in unc_sic_files:
            with np.load(unc_sic_file) as ds:
                unc_myi.append(ds['unc_sic'])
        return unc_myi

    def get_unc_sid_files(self, idate):
        unc_sid_files = []
        for i in range(7):
            src_year = idate.year - i
            unc_sid_file = f'{self.unc_dir}/{src_year}/unc_sid_{src_year}0905_{idate.year}{idate.month:02d}{idate.day:02d}.npz'
            if os.path.exists(unc_sid_file):
                unc_sid_files.append(unc_sid_file)
        return unc_sid_files

    def load_unc_sid(self, idate, unc_obs):
        unc_sid_files = self.get_unc_sid_files(idate)
        unc_sid = [np.zeros_like(unc_obs)]
        for unc_sid_file in unc_sid_files:
            with np.load(unc_sid_file) as ds:
                unc_sid.append(ds['unc_sid'])
        return unc_sid
    
    def get_unc_sidsic_files(self, idate):
        unc_sid_files = []
        for i in range(7):
            src_year = idate.year - i
            unc_sid_file = f'{self.unc_dir}/{src_year}/unc_sidsic_{src_year}0905_{idate.year}{idate.month:02d}{idate.day:02d}.npz'
            if os.path.exists(unc_sid_file):
                unc_sid_files.append(unc_sid_file)
        return unc_sid_files
    
    def load_unc_sidsic(self, idate, unc_obs):
        unc_sid_files = self.get_unc_sidsic_files(idate)
        unc_sidsic = [np.zeros_like(unc_obs)]
        for unc_sid_file in unc_sid_files:
            with np.load(unc_sid_file) as ds:
                unc_sidsic.append(ds['unc_sic_sid'])
        return unc_sidsic

    def compute_unc_frac(self, unc_obs, unc_myi):
        unc_fyi = np.hypot(unc_obs, unc_myi[0])
        unc_frac = [unc_fyi]
        for i in range(len(unc_myi)):
            unc_square_sum = 0
            for j in range(min(i+2, len(unc_myi))):
                unc_square_sum += unc_myi[j]**2
            unc_frac.append(np.sqrt(unc_square_sum))
        return unc_frac

    def compute_unc_tot(self, unc_frac, unc_sic_sid):
        unc_tot = []
        unc_age = 0
        for i in range(len(unc_frac)):
            unc_tot.append(np.hypot(unc_frac[i], unc_sic_sid[i]))
            unc_age += (i * unc_tot[i]/100)**2
        unc_age = np.sqrt(unc_age)
        return unc_tot, unc_age

    def compute_uncertainty(self, idate):
        unc_tot_file = f'{self.unc_dir}/{idate.year}/unc_tot_{idate.year}{idate.month:02d}{idate.day:02d}.npz'
        if os.path.exists(unc_tot_file):
            return
        unc_obs = self.load_unc_obs(idate)
        unc_myi = self.load_unc_myi(idate)
        unc_sidsic = self.load_unc_sidsic(idate, unc_obs)
        unc_frac = self.compute_unc_frac(unc_obs, unc_myi)
        unc_tot, unc_age = self.compute_unc_tot(unc_frac, unc_sidsic)
        np.savez(unc_tot_file, unc_tot=np.array(unc_tot).astype(np.float16), unc_age=unc_age.astype(np.float16))
