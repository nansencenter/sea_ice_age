from datetime import datetime, timedelta
import glob
import os

from cartopy.crs import NorthPolarStereo, LambertAzimuthalEqualArea
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt, gaussian_filter
from tqdm import tqdm

class SID_loader_CDR:
    def load_data(self, sid_files):
        ddx = []
        ddy = []
        unc = []
        for sid_file in tqdm(sid_files):
            with Dataset(sid_file) as dds:
                dx = dds['dX'][0].filled(np.nan)
                dy = dds['dY'][0].filled(np.nan)
                un = dds['uncert_dX_and_dY'][0].filled(np.nan)
            ddx.append(dx)
            ddy.append(dy)
            unc.append(un)
        ddx = np.array(ddx)
        ddy = np.array(ddy)
        unc = np.array(unc)
        return ddx, ddy, unc

def get_si_conc_for_sid_dates(sid_dates, sic_data_dir):
    sic_file_format = f'{sic_data_dir}/%Y/%m/ice_conc_nh_ease2-*_%Y%m%d1200.nc'
    sic_files = [date.strftime(sic_file_format) for date in sid_dates]
    print(len(sic_files), sic_files[0], sic_files[-1])

    sic_all = []
    unc_all = []
    valid_sid_idx = []
    for i, sic_file in tqdm(enumerate(sic_files), total=len(sid_dates)):
        match_files = glob.glob(sic_file)
        if len(match_files) > 0:
            with Dataset(match_files[0]) as dds:
                conc = dds['ice_conc'][0].filled(np.nan)
                unct = dds['algorithm_standard_uncertainty'][0].filled(np.nan)
                sic_all.append(conc)
                unc_all.append(unct)
                valid_sid_idx.append(i)
    sic = np.array(sic_all)
    sic[sic < 0] = np.nan
    unc = np.array(unc_all)
    unc[unc < 0] = np.nan
    return sic, unc, valid_sid_idx

def filter_3d_array(array, ice_mask, distance=5, gauss_filter_kernel=(0.5, 1., 1.), truncate=2):
    mask = np.isnan(array)
    dist, indi = distance_transform_edt(mask, return_distances=True, return_indices=True)
    gpi = dist <= distance
    r,c,i = indi[:, gpi]
    array_f1 = np.array(array)
    array_f1[gpi] = array[r,c,i]
    array_f2 = gaussian_filter(array_f1, gauss_filter_kernel, truncate=truncate)
    array_f2[ice_mask] = np.nan
    return array_f2

def save_data(out_data_dir, sid_dates, ddxf00, ddyf00, sid_unc_f00, si_con, sic_unc):
    for i, date in tqdm(enumerate(sid_dates), total=len(sid_dates)):
        osubdir = date.strftime(f'{out_data_dir}/%Y')
        os.makedirs(osubdir, exist_ok=True)
        ofile = f'{osubdir}/ice_drift_nh_ease2-750_cdr-v1p0_24h-{date.strftime("%Y%m%d")}1200.nc.npz'
        np.savez(ofile, u=ddxf00[i], v=-ddyf00[i], sid_unc=sid_unc_f00[i], c=si_con[i], sic_unc=sic_unc[i])


class SID_loader_iCDR:
    srs_src = NorthPolarStereo(-45, 70)
    srs_dst = LambertAzimuthalEqualArea(0, 90)
    def __init__(self, sid_dir_src=None, sid_dir_dst=None):
        self.sid_file_src = f'{sid_dir_src}/2020/12/ice_drift_nh_polstere-625_multi-oi_202012291200-202012311200.nc'
        self.sid_file_dst = f'{sid_dir_dst}/2020/12/ice_drift_nh_ease2-750_cdr-v1p0_24h-202012301200.nc'

    def get_transformation_grids(self):
        with Dataset(self.sid_file_src) as dds:
            #dx_src = dds['dX'][0]
            #dy_src = dds['dY'][0]
            sid_x_src = dds['xc'][:] * 1000
            sid_y_src = dds['yc'][:] * 1000

        with Dataset(self.sid_file_dst) as dds:
            #dx_dst = dds['dX'][0]
            #dy_dst = dds['dY'][0]
            sid_x_dst = dds['xc'][:] * 1000
            sid_y_dst = dds['yc'][:] * 1000

        sid_x_dst_grd, sid_y_dst_grd = np.meshgrid(sid_x_dst, sid_y_dst)
        tmp = self.srs_src.transform_points(self.srs_dst, sid_x_dst_grd, sid_y_dst_grd)
        sid_x_dst_grd_sid, sid_y_dst_grd_sid = tmp[:,:,0], tmp[:,:,1]

        # rotate vectors
        am = np.array([
            [np.cos(np.radians(-45)), -np.sin(np.radians(-45))],
            [np.sin(np.radians(-45)),  np.cos(np.radians(-45))],
        ])
        self.sid_x_src = sid_x_src
        self.sid_y_src = sid_y_src
        self.sid_x_dst_grd_sid = sid_x_dst_grd_sid
        self.sid_y_dst_grd_sid = sid_y_dst_grd_sid
        self.am = am

    def get_sid405_on_laea(self, sid_file):
        # READ ALL SID DATA (WITH ROTATION and REPROJECTION)
        with Dataset(sid_file) as dds:
            dx_src = dds['dX'][0]
            dy_src = dds['dY'][0]
            uncert = dds['uncert_dX_and_dY'][0].filled(0)
        dx_src_rot = (dx_src * self.am[0,0] + dy_src * self.am[0,1]) / 2
        dy_src_rot = (dx_src * self.am[1,0] + dy_src * self.am[1,1]) / 2
        rgi_dx = RegularGridInterpolator(
            (self.sid_y_src[::-1], self.sid_x_src), dx_src_rot[::-1], method='nearest', bounds_error=False)
        rgi_dy = RegularGridInterpolator(
            (self.sid_y_src[::-1], self.sid_x_src), dy_src_rot[::-1], method='nearest', bounds_error=False)
        rgi_un = RegularGridInterpolator(
            (self.sid_y_src[::-1], self.sid_x_src), uncert[::-1], method='nearest', bounds_error=False)
        dx_src_dst = rgi_dx((self.sid_y_dst_grd_sid, self.sid_x_dst_grd_sid))
        dy_src_dst = rgi_dy((self.sid_y_dst_grd_sid, self.sid_x_dst_grd_sid))
        uncert_dst = rgi_un((self.sid_y_dst_grd_sid, self.sid_x_dst_grd_sid))
        return dx_src_dst, dy_src_dst, uncert_dst
    
    def load_data(self, sid_files):
        self.get_transformation_grids()
        ddx = []
        ddy = []
        unc = []
        for sid_file in tqdm(sid_files):
            dx, dy, un = self.get_sid405_on_laea(sid_file)
            ddx.append(dx.filled(np.nan))
            ddy.append(dy.filled(np.nan))
            unc.append(un)
        ddx = np.array(ddx)
        ddy = np.array(ddy)
        unc = np.array(unc)
        return ddx, ddy, unc

def load_and_save(year, sid_dir, sic_data_dir, out_data_dir, loader, min_sic=15, distance=10, dst_date_offset=0):
    sid_files = sorted(glob.glob(f'{sid_dir}/{year}/??/*.nc'))
    sid_dates = [f.split('-')[-1][:8] for f in sid_files]
    sid_dates = [datetime.strptime(d, '%Y%m%d') + timedelta(dst_date_offset) for d in sid_dates]
    print(len(sid_dates), sid_dates[0], sid_dates[-1])
    ddx, ddy, sid_unc = loader.load_data(sid_files)
    si_con, sic_unc, valid_sid_idx = get_si_conc_for_sid_dates(sid_dates, sic_data_dir=sic_data_dir)
    if len(valid_sid_idx) == 0:
        print(f'No valid SIC data for {year}')
        return
    ice_mask = (si_con[:, 1:-1:3, 1:-1:3] < min_sic) + np.isnan(si_con[:, 1:-1:3, 1:-1:3])
    ddxf00 = filter_3d_array(ddx[valid_sid_idx], ice_mask, distance=distance)
    ddyf00 = filter_3d_array(ddy[valid_sid_idx], ice_mask, distance=distance)
    sid_unc_f00 = filter_3d_array(sid_unc[valid_sid_idx], ice_mask, distance=distance)
    valid_sid_dates = [sid_dates[i] for i in valid_sid_idx]
    save_data(out_data_dir, valid_sid_dates, ddxf00, ddyf00, sid_unc_f00, si_con, sic_unc)