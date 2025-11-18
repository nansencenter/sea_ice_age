import glob
from datetime import datetime, timedelta
import os

import numpy as np
from scipy.interpolate import RectBivariateSpline
from skimage.util import view_as_windows
from scipy.ndimage import uniform_filter

from lmsiage.mesh_file import MeshFile
from lmsiage.utils import fill_gaps
from lmsiage.zarr_index import get_file_arrays


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

def compute_std(ifile, col, row):
    myi = np.load(ifile)['c']
    myi_grd = np.zeros((row.max() + 1, col.max() + 1)) + np.nan
    myi_grd[row, col] = myi
    std_grd = {}
    std_sizes = list(range(3,31,3))
    for size in std_sizes:
        std_grd[size] = std_filter_using_windows(myi_grd, size=size).astype(np.float16)
        std_grd[size][np.isnan(std_grd[size])] = 0
    return std_grd


class ComputeSidSicUncertainty:
    def __init__(self, sia_dir):
        self.sia_dir = sia_dir

    def save_sid_sic_uncertainty(self, idate):
        xel_min = -2000
        yel_min = -2450
        node_distance_mean = 30

        afile = f'{self.sia_dir}/age/{idate.year}/age_{idate.strftime("%Y%m%d")}.npz'
        if not os.path.exists(afile):
            print(f'File {afile} does not exist.')
            return
        with np.load(afile) as ds:
            x = ds['x']
            y = ds['y']
            t = ds['t'].astype(int)
        xel = x[t].mean(axis=1)
        yel = y[t].mean(axis=1)
        col = ((xel  - xel_min) / 25).astype(int)
        row = ((yel  - yel_min) / 25).astype(int)

        ifiles = []
        unc_sid_files = []
        for i in range(7):
            src_year = idate.year - i
            ifile = f'{self.sia_dir}/sic/{src_year}/sic_{src_year}0915_{idate.strftime("%Y%m%d")}.npz'
            unc_sid_file = f'{self.sia_dir}/unc/{src_year}/unc_sid_{src_year}0905_{idate.strftime("%Y%m%d")}.npz'

            if os.path.exists(ifile) and os.path.exists(unc_sid_file):
                ifiles.append(ifile)
                unc_sid_files.append(unc_sid_file)

        for ifile, unc_sid_file in zip(ifiles, unc_sid_files):
            ofile = unc_sid_file.replace('unc_sid_', 'unc_sidsic_')
            if os.path.exists(ofile):
                continue
            std_grd = compute_std(ifile, col, row)
            unc_sid_sum = np.load(unc_sid_file)['unc_sid']
            min_rel_size_factor = 3
            min_rel_sizes = np.arange(3, 31, 3)
            sid_unc_rel_to_el_size = unc_sid_sum / node_distance_mean
            unc_sic_sid = np.zeros_like(unc_sid_sum)
            for min_rel_size in min_rel_sizes:
                uncert_elems = np.nonzero(sid_unc_rel_to_el_size > min_rel_size)[0]
                if uncert_elems.size == 0:
                    continue
                col_unc_sid = col[uncert_elems]
                row_unc_sid = row[uncert_elems]
                unc_sic_sid[uncert_elems] = std_grd[min_rel_size][row_unc_sid, col_unc_sid]
            np.savez(ofile, unc_sic_sid=unc_sic_sid.astype(np.float16))


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
            ar = d['ar']
            ar[ar == 0] = 0.01
        return sic_unc, x, y, t, src2dst, weights, ar

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

            mesh_file, _ = get_mesh_files(idate, self.mesh_dir, self.mesh_init_file)
            unc_sic_grd, x, y, t, src2dst, weights, ar = self.load_sic_data(ifile, mesh_file)
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
                unc_sic_sum_pro /= ar
                unc_sic_sum_pro = np.clip(unc_sic_sum_pro, 0, 100)
                
                # advect sic_min from previous time step
                sic_min_pro = np.zeros(src2dst[:,1].max()+1)
                np.add.at(sic_min_pro, src2dst[:,1], sic_min[src2dst[:,0]] * weights)
                sic_min_pro /= ar
                sic_min_pro = np.clip(sic_min_pro, 0, 100)

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
    def __init__(self, mesh_dir, unc_dir, ifiles, idates, n_steps, size, xc, yc):
        self.mesh_dir = mesh_dir
        self.unc_dir = unc_dir
        self.ifiles = ifiles
        self.idates = idates
        self.n_steps = n_steps
        self.size = size
        self.xc = xc
        self.yc = yc

    def load_sid_data(self, ifile, mesh_file):
        with np.load(ifile) as d:
            sid_unc = d['sid_unc']
        
        mf = MeshFile(mesh_file)
        x, y, t = mf.load(['x', 'y', 't'], as_dict=False)
        return sid_unc, x, y, t
    
    def compute_smearing_uncertainty(self, year):
        """ Compute smearing uncertainty and integrated drift uncertainty """
        unc_sid_name = 'unc_sid'
        mesh_files = sorted(glob.glob(f'{self.mesh_dir}/{year}/mesh_*.zip'))
        mesh_dates = [datetime.strptime(os.path.basename(mesh_file), 'mesh_%Y%m%d.zip') + timedelta(days=0.5) for mesh_file in mesh_files]
        unc_files = [mesh_date.strftime(f'{self.unc_dir}/%Y/unc_%Y%m%d.zip') for mesh_date in mesh_dates]
            
        dst_dir = f'{self.unc_dir}/{year}'
        os.makedirs(dst_dir, exist_ok=True)

        for mesh_file, mesh_date, unc_file in zip(mesh_files, mesh_dates, unc_files):
            if unc_sid_name in get_file_arrays(unc_file):
                continue

            ifile = self.ifiles[self.idates.index(mesh_date)]
            unc_sid_grd, x, y, t = self.load_sid_data(ifile, mesh_file)
            unc_sid_fil = fill_gaps(unc_sid_grd, np.isnan(unc_sid_grd), distance=144)
            # compute uncertainty of drift field after smoothing
            unc_sid_smt = uniform_filter(unc_sid_fil**2, size=self.size) / self.size
            unc_sid_smt[np.isnan(unc_sid_grd)] = 0
            # interpolate uncertainty to mesh
            unc_sid = RectBivariateSpline(self.xc[1:-1:3], self.yc[1:-1:3], unc_sid_smt[::-1], kx=1, ky=1)(y[t].mean(axis=1), x[t].mean(axis=1), grid=False)
            # Save uncertainty
            unc_mf = MeshFile(unc_file)
            unc_mf.save({unc_sid_name: unc_sid.astype(np.float16)}, mode='a')

    def compute_integrated_uncertainty(self, unc_sid_files, unc_sid_dates, start_idx):
        """ Compute integrated drift uncertainty """
        stop_idx = min(start_idx + self.n_steps, len(unc_sid_files))
        unc_sid_int = None
        start_date = unc_sid_dates[start_idx]
        unc_sid_name = f'unc_sid{start_date.year}'
        print('Processing', start_date, unc_sid_dates[stop_idx-1])
        for i in range(start_idx, stop_idx):
            unc_date = unc_sid_dates[i]
            mesh_file = unc_date.strftime(f'{self.mesh_dir}/%Y/mesh_%Y%m%d.zip')
            unc_file = unc_date.strftime(f'{self.unc_dir}/%Y/unc_%Y%m%d.zip')

            mf = MeshFile(unc_file)
            # load precomputed integared uncertainty, skip the rest of processing
            if unc_sid_name in get_file_arrays(unc_file):
                unc_sid_int = mf.load([unc_sid_name], as_dict=False)[0]
                continue

            # load smearing uncertainty
            unc_sid = mf.load(['unc_sid'], as_dict=False)[0]
            if unc_sid_int is None:
                # first time step: initialize integrated uncertainty
                unc_sid_int = np.array(unc_sid)
            else:
                # later steps: advect unc_sid_int from previous time step
                src2dst, weights = MeshFile(mesh_file).load(['src2dst', 'weights'], as_dict=False)
                unc_sid_int_advected = np.zeros(src2dst[:,1].max()+1)
                np.add.at(unc_sid_int_advected, src2dst[:,1], unc_sid_int[src2dst[:,0]] * weights)
                # compute new integrated uncertainty after advection
                unc_sid_int = np.hypot(unc_sid, unc_sid_int_advected)
            
            # Save integrated uncertainty
            MeshFile(unc_file).save({unc_sid_name: unc_sid_int.astype(np.float16)}, mode='a')
        print('Done', start_date, self.idates[stop_idx-1])



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
        return unc_obs, len(f)

    def get_unc_sic_files(self, idate):
        unc_sic_files = []
        for i in range(7):
            if i == 0 and idate.month == 9 and idate.day < 15:
                continue
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

    def compute_unc_frac(self, unc_obs, unc_myi, len_f):
        unc_fyi = np.hypot(unc_obs, unc_myi[0])
        unc_frac = [unc_fyi]
        for i in range(len_f-1):
            unc_square_sum = 0
            for j in range(min(i+2, len(unc_myi))):
                unc_square_sum += unc_myi[j]**2
            unc_frac.append(np.sqrt(unc_square_sum))
        return unc_frac

    def compute_unc_tot(self, unc_frac, unc_sic_sid):
        unc_tot = []
        unc_age = 0
        for i in range(len(unc_frac)):
            try:
                unc_tot.append(np.hypot(unc_frac[i], unc_sic_sid[i]))
            except IndexError:
                return None, None
            unc_age += (i * unc_tot[i]/100)**2
        unc_age = np.sqrt(unc_age)
        return unc_tot, unc_age

    def compute_uncertainty(self, idate):
        unc_tot_file = f'{self.unc_dir}/{idate.year}/unc_tot_{idate.year}{idate.month:02d}{idate.day:02d}.npz'
        if os.path.exists(unc_tot_file):
            return
        unc_obs, len_f = self.load_unc_obs(idate)
        unc_myi = self.load_unc_myi(idate)
        unc_sidsic = self.load_unc_sidsic(idate, unc_obs)
        unc_frac = self.compute_unc_frac(unc_obs, unc_myi, len_f)
        unc_tot, unc_age = self.compute_unc_tot(unc_frac, unc_sidsic)
        if unc_tot is None:
            print(f'Uncertainty for {idate.strftime("%Y%m%d")} is None')
            raise
        np.savez(unc_tot_file, unc_tot=np.array(unc_tot).astype(np.float16), unc_age=unc_age.astype(np.float16))
