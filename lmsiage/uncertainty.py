import glob
from datetime import datetime, timedelta
import os
import warnings

import numpy as np
from scipy.interpolate import RectBivariateSpline
from skimage.util import view_as_windows
from scipy.ndimage import uniform_filter
from tqdm import tqdm, trange

from lmsiage.mesh_file import MeshFile
from lmsiage.utils import fill_gaps
from lmsiage.zarr_index import get_file_arrays

warnings.filterwarnings('ignore')

class ComputeObsUncertainty:
    """ Compute observation uncertainties of SIC and SID mapped to mesh grid."""
    def __init__(self, obs_files, obs_dates, mesh_dir, unc_dir, xc, yc, size=3):
        self.obs_files = obs_files
        self.obs_dates = obs_dates
        self.mesh_dir = mesh_dir
        self.unc_dir = unc_dir
        self.xc = xc
        self.yc = yc
        self.size = size

    def load_data(self, obs_file, mesh_file):
        """ Load SIC and SID uncertainty from observation file and mesh grid info from mesh file."""
        with np.load(obs_file) as d:
            sic_unc = d['sic_unc']
            sid_unc = d['sid_unc']
        
        mf = MeshFile(mesh_file)
        x, y, t = mf.load(['x', 'y', 't'], as_dict=False)
        return sic_unc, sid_unc, x, y, t

    def compute_obs_uncertainty(self, year):
        """ Compute observation uncertainties of SIC and SID mapped to mesh grid for a given year."""
        mesh_files = sorted(glob.glob(f'{self.mesh_dir}/{year}/mesh_*.zip'))
        mesh_dates = [datetime.strptime(os.path.basename(mesh_file), 'mesh_%Y%m%d.zip') + timedelta(days=0.5) for mesh_file in mesh_files]
        unc_files = [mesh_date.strftime(f'{self.unc_dir}/%Y/unc_%Y%m%d.zip') for mesh_date in mesh_dates]
        print(f'Computing OBS uncertainty for {len(mesh_files)} mesh files for year {year}')
        dst_dir = f'{self.unc_dir}/{year}'
        os.makedirs(dst_dir, exist_ok=True)

        unc_sic_name = 'unc_sic'
        unc_sid_name = 'unc_sid'

        for mesh_file, mesh_date, unc_file in zip(tqdm(mesh_files), mesh_dates, unc_files):
            file_arrays = get_file_arrays(unc_file)
            if unc_sic_name in file_arrays and unc_sid_name in file_arrays:
                continue
            
            obs_file = self.obs_files[self.obs_dates.index(mesh_date)]
            # load SIC and SID uncert from obs file
            unc_sic_grd, unc_sid_grd, x, y, t = self.load_data(obs_file, mesh_file)

            # process SIC uncertainty
            unc_sic_fil = fill_gaps(unc_sic_grd, np.isnan(unc_sic_grd), distance=144)
            unc_sic = RectBivariateSpline(self.xc, self.yc, unc_sic_fil[::-1], kx=1, ky=1)(y[t].mean(axis=1), x[t].mean(axis=1), grid=False)

            # process SID uncertainty to mesh
            unc_sid_fil = fill_gaps(unc_sid_grd, np.isnan(unc_sid_grd), distance=144)
            # compute uncertainty of drift field after smoothing
            unc_sid_smt = uniform_filter(unc_sid_fil**2, size=self.size) / self.size
            unc_sid_smt[np.isnan(unc_sid_grd)] = 0
            unc_sid = RectBivariateSpline(self.xc[1:-1:3], self.yc[1:-1:3], unc_sid_smt[::-1], kx=1, ky=1)(y[t].mean(axis=1), x[t].mean(axis=1), grid=False)

            # Save uncertainty
            unc_mf = MeshFile(unc_file)
            unc_mf.save({
                unc_sic_name: unc_sic.astype(np.float16),
                unc_sid_name: unc_sid.astype(np.float16),
                }, mode='a')


class ComputeIntegratedUncertainty:
    def __init__(self, mesh_dir, age_dir, unc_dir, n_steps):
        self.n_steps = n_steps
        self.mesh_dir = mesh_dir
        self.age_dir = age_dir
        self.unc_dir = unc_dir

    def compute(self, unc_sic_dates, start_idx):
        """ Compute integrated uncertainties over n_steps starting from start_idx."""
        # init integrated uncerts
        sic_min = None
        unc_sic_min = None
        unc_sid_int = None

        stop_idx = min(start_idx + self.n_steps, len(unc_sic_dates))
        start_date = unc_sic_dates[start_idx]
        
        unc_sic_name = f'unc_sic{start_date.year}'
        sic_min_name = f'sic_min{start_date.year}'
        sic_name = f'sic{start_date.strftime("%Y%m%d")}'
        unc_sid_name = f'unc_sid{start_date.year}'

        print(f'Computing uncertainty for {start_date.strftime("%Y%m%d")}, {unc_sic_dates[stop_idx-1]}')
        for i in trange(start_idx, stop_idx):
            unc_date = unc_sic_dates[i]
            mesh_file = unc_date.strftime(f'{self.mesh_dir}/%Y/mesh_%Y%m%d.zip')
            age_file = unc_date.strftime(f'{self.age_dir}/%Y/age_%Y%m%d.zip')
            unc_file = unc_date.strftime(f'{self.unc_dir}/%Y/unc_%Y%m%d.zip')
            
            file_arrays = get_file_arrays(unc_file)
            if unc_sic_name in file_arrays and unc_sid_name in file_arrays:
                # load precomputed min concentration and its uncertainty, skip the rest of processing
                sic_min, unc_sic_min, unc_sid_int = MeshFile(unc_file).load([sic_min_name, unc_sic_name, unc_sid_name], as_dict=False)
                continue
            
            age_arrays = get_file_arrays(age_file)
            if sic_name not in age_arrays:
                # no sic data for this date, stop processing
                break

            # load obs concentration, uncertainty
            sic = MeshFile(age_file).load([sic_name], as_dict=False)[0]
            unc_sic, unc_sid = MeshFile(unc_file).load(['unc_sic', 'unc_sid'], as_dict=False)
            if unc_sic_min is None:
                # first time step: initialize minimal concentration and its uncertainty
                unc_sic_min = np.array(unc_sic)
                sic_min = np.array(sic)
                # initialize drift integrated uncertainty
                unc_sid_int = np.array(unc_sid)
            else:
                # later steps: advect sic_min, unc_sic_min from previous time step
                src2dst, weights = MeshFile(mesh_file).load(['src2dst', 'weights'], as_dict=False)

                unc_sic_min_advected = np.zeros(src2dst[:,1].max()+1)
                np.add.at(unc_sic_min_advected, src2dst[:,1], unc_sic_min[src2dst[:,0]] * weights)
                sic_min_advected = np.zeros(src2dst[:,1].max()+1)
                np.add.at(sic_min_advected, src2dst[:,1], sic_min[src2dst[:,0]] * weights)
                unc_sid_int_advected = np.zeros(src2dst[:,1].max()+1)
                np.add.at(unc_sid_int_advected, src2dst[:,1], unc_sid_int[src2dst[:,0]] * weights)

                # keep minimal value of SIC
                sic_min = np.min([sic_min_advected, sic], axis=0)

                # keep uncertainty of the minimal value of SIC
                min_sic_ids = sic < sic_min_advected
                unc_sic_min = unc_sic_min_advected
                unc_sic_min[min_sic_ids] = unc_sic[min_sic_ids]

                # compute new drift integrated uncertainty after advection
                unc_sid_int = np.hypot(unc_sid, unc_sid_int_advected)
            
            # Save minimal concentration and its uncertainty
            MeshFile(unc_file).save({
                    sic_min_name: sic_min.astype(np.float16),
                    unc_sic_name: unc_sic_min.astype(np.float16),
                    unc_sid_name: unc_sid_int.astype(np.float16),
                 }, mode='a')
        print('Done', start_date, unc_sic_dates[stop_idx-1])


class ComputeAgeUncertainty:
    MAX_LEN_FRACTIONS = 6

    """ Compute sea ice age uncertainty from sea ice concentration and sea ice drift uncertainties."""
    def __init__(self, mesh_dir, age_dir, unc_dir, save_names=('unc_age', 'unc_fracs'), force=False):
        self.mesh_dir = mesh_dir
        self.age_dir = age_dir
        self.unc_dir = unc_dir
        self.save_names = save_names
        self.force = force

    def __call__(self, unc_file):
        """ Process a single uncertainty file to compute age uncertainty."""
        mesh_file, age_file, unc_file = self.get_files(unc_file)
        status = self.load_data(mesh_file, age_file, unc_file)
        if status < 0:
            return status
        self.compute()
        self.save(unc_file)
        return 0
    
    def get_files(self, unc_file):
        """ Get corresponding mesh and age files for a given uncertainty file."""
        date = datetime.strptime(os.path.basename(unc_file), 'unc_%Y%m%d.zip')
        mesh_file = date.strftime(f'{self.mesh_dir}/%Y/mesh_%Y%m%d.zip')
        age_file = date.strftime(f'{self.age_dir}/%Y/age_%Y%m%d.zip')
        return mesh_file, age_file, unc_file

    def load_data(self, mesh_file, age_file, unc_file, xel_min=-2000, yel_min=-2450):
        """ Load mesh, age, and uncertainty data from files."""
        x, y, t = MeshFile(mesh_file).load(['x', 'y', 't'], as_dict=False)
        xel = x[t].mean(axis=1)
        yel = y[t].mean(axis=1)
        self.col = ((xel  - xel_min) / 25).astype(int)
        self.row = ((yel  - yel_min) / 25).astype(int)

        unc_names = get_file_arrays(unc_file)
        self.unc_years = sorted([n[7:] for n in unc_names if 'unc_sid' in n and len(n) > 7])[-self.MAX_LEN_FRACTIONS:]
        if len(self.unc_years) == 0:
            # no uncertainty years found
            return -1
        
        save_names_in_file = all(name in unc_names for name in self.save_names)
        if save_names_in_file and not self.force:
            # all data already computed
            return -2
        
        self.sic_myi = {}
        self.unc_myi = {}
        self.unc_sid = {}
        for uy in self.unc_years:
            self.sic_myi[uy], self.fracs = MeshFile(age_file).load([f'sic{uy}0915', 'f'], as_dict=False)
            self.unc_sic, self.unc_myi[uy], self.unc_sid[uy] = MeshFile(unc_file).load(['unc_sic', f'unc_sic{uy}', f'unc_sid{uy}'], as_dict=False)
        return 0

    def compute(self):
        """ Compute all uncertainties."""
        self.std_grd = {}
        self.unc_sic_sid = {}
        self.unc_com = {}
        for uy in self.unc_years:
            self.std_grd[uy] = self.compute_std_grids(uy)
            self.unc_sic_sid[uy] = self.compute_sid_sic_uncert(uy)
            self.unc_com[uy] = np.sqrt(self.unc_myi[uy]**2 + self.unc_sic_sid[uy]**2)
        self.unc_fracs = self.compute_fraction_uncertainties()
        self.unc_age = self.compute_age_uncertainty()

    def compute_std_grids(self, uy, min_size=3, max_size=31, step=3):
        """ Compute standard deviation grids for multiple window sizes."""
        myi_grd = np.zeros((self.row.max() + 1, self.col.max() + 1)) + np.nan
        myi_grd[self.row, self.col] = self.sic_myi[uy]
        std_grd = {}
        std_sizes = list(range(min_size, max_size, step))
        for size in std_sizes:
            std_grd[size] = self.std_filter_using_windows(myi_grd, size=size).astype(np.float16)
            std_grd[size][np.isnan(std_grd[size])] = 0
        return std_grd

    def std_filter_using_windows(self, image, size=25):
        """ Compute the standard deviation of an image using sliding windows. """
        pad_width = size // 2  # for 5x5 window
        image_pad = np.pad(image, pad_width, mode='constant', constant_values=np.nan)
        # Create a sliding window view of the image
        windows = view_as_windows(image_pad, (size, size), step=1)
        # Compute the standard deviation for each window
        std_image = np.nanstd(windows, axis=(2, 3))
        return std_image

    def compute_sid_sic_uncert(self, uy, node_distance_mean=30, min_size=3, max_size=31, step=3):
        """ Select standard deviation values based on relative element sizes."""
        min_rel_sizes = np.arange(min_size, max_size, step)
        sid_unc_rel_to_el_size = self.unc_sid[uy] / node_distance_mean
        unc_sic_sid = np.zeros_like(self.unc_sid[uy])
        for min_rel_size in min_rel_sizes:
            uncert_elems = np.nonzero(sid_unc_rel_to_el_size > min_rel_size)[0]
            if uncert_elems.size == 0:
                continue
            col_unc_sid = self.col[uncert_elems]
            row_unc_sid = self.row[uncert_elems]
            unc_sic_sid[uncert_elems] = self.std_grd[uy][min_rel_size][row_unc_sid, col_unc_sid]
        return unc_sic_sid

    def compute_fraction_uncertainties(self):
        """ Compute fractional uncertainties for all source year combinations."""
        # uncerts of advected MYI as a list sorted by year
        # min of advected
        unc_myi = [self.unc_myi[key] for key in sorted(self.unc_myi.keys())]
        # combined
        unc_com = [self.unc_com[key] for key in sorted(self.unc_com.keys())]
        # uncerts of ice age fractions
        # sigma_NYI^2 = sigma_COM,N-2^2 + sigma_COM,N-1^2
        # sigma_FYI^2 = sigma_MYI,N-1^2 + sigma_SIC^2
        unc_fracs = [unc_com[0]] + [np.hypot(unc_com[i], unc_com[i+1]) for i in range(len(unc_com) - 1)] + [np.hypot(unc_myi[-1], self.unc_sic)]
        return np.array(unc_fracs)

    def compute_age_uncertainty(self):
        """ Compute age uncertainties from fractional uncertainties and fractions.
        (sigma_q / q)^2 = (sigma_x / x^2) + (sigma_y / y^2)
        where x = sum(ages_vec * fracs), y = sum(fracs)
        """
        ages_vec = (1 + np.arange(len(self.unc_fracs), dtype=np.float32)[::-1])[:, None]

        sigma_x = np.sum(ages_vec * self.unc_fracs ** 2, axis=0)
        sigma_y = np.sum(self.unc_fracs ** 2, axis=0)
        x_ = np.sum(ages_vec * self.fracs, axis=0)
        y_ = np.sum(self.fracs, axis=0)
        x_[y_ < 1] = np.nan
        y_[y_ < 1] = np.nan
        unc_age = np.sqrt(sigma_x / x_**2 + sigma_y / y_**2) * x_ / y_
        unc_age[np.isnan(unc_age)] = 0
        return unc_age

    def save(self, unc_file):
        """ Save computed uncertainties to file."""
        data = {}
        for i in self.save_names:
            if isinstance(self.__dict__[i], dict):
                save_items = {f'{i}{k}': self.__dict__[i][k] for k in self.__dict__[i]}
            else:
                save_items = {i: self.__dict__[i] }
            data.update(save_items)
        mf = MeshFile(unc_file)
        mf.save(data, mode='o')
