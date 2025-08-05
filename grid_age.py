from datetime import datetime, timedelta
import os

from netCDF4 import Dataset
from geodataset.geodataset import GeoDatasetWrite
import numpy as np
import pyproj

from utils import IrregularGridInterpolator

class GridAge:
    def __init__(self, xgrd, ygrd, mask, force=False):
        self.xgrd = xgrd
        self.ygrd = ygrd
        self.mask = mask
        self.force = force

    def __call__(self, age_file):
        age_grd_file = age_file.replace('.npz', '_grd.npz').replace('/age/', '/age_grd/')
        unc_file = age_file.replace('/age/', '/unc/').replace('/age_', '/unc_tot_')
        if os.path.exists(age_grd_file) and not self.force:
            return
        try:
            with np.load(age_file) as data:
                t = data['t']
                x = data['x']
                y = data['y']
                a = data['a']
                f = data['f']
        except:
            raise ValueError(f'Cannot load from {age_file}')
        
        try:
            with np.load(unc_file) as data:                
                unc_tot = data['unc_tot']
                unc_age = data['unc_age']
        except:
            raise ValueError(f'Cannot load from {unc_file}')

        try:
            igi = IrregularGridInterpolator(x, y, self.xgrd, self.ygrd, t)
        except:
            raise ValueError(f'Cannot create IGI {age_file}')

        src_data = np.vstack([a[None], f, unc_age[None], unc_tot[:len(f)]])
        dst_data = []
        for d in src_data:
            try:
                dgrd = igi.interp_field(d)
            except:
                import ipdb; ipdb.set_trace()
                raise ValueError(f'Cannot interpolate {age_file}')
            dgrd[self.mask == 0] = np.nan
            dst_data.append(dgrd.astype(np.float32))
        
        frac_num = len(f)
        save_data_names = (
            ['age'] + 
            [f'fraction_{frac_num - i}yi' for i in range(frac_num)] +
            ['age_unc'] + 
            [f'fraction_{i}yi_unc' for i in range(1, len(f) + 1)]
        )
        if len(dst_data) != len(save_data_names):
            raise ValueError(f'Length mismatch: {len(dst_data)} != {len(save_data_names)} for {age_file}')
        save_data = {}
        for i, name in enumerate(save_data_names):
            save_data[name] = dst_data[i]
        dst_date_dir = os.path.split(age_grd_file)[0]
        os.makedirs(dst_date_dir, exist_ok=True)
        np.savez_compressed(age_grd_file, **save_data)

class SeaIceAgeDataset(GeoDatasetWrite):
    """ wrapper for netCDF4.Dataset with info about Ice Age products """
    grid_mapping_variable = 'Lambert_Azimuthal_Equal_Area'
    projection = pyproj.Proj("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
    global_attributes_source = None
    global_attributes_title = None

    def get_grid_mapping_ncattrs(self):
        return dict(
            grid_mapping_name = "lambert_azimuthal_equal_area" ,
            longitude_of_projection_origin = 0. ,
            latitude_of_projection_origin = 90. ,
            false_easting = 0. ,
            false_northing = 0. ,
            semi_major_axis = 6378137. ,
            inverse_flattening = 298.257223563 ,
            proj4_string = "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs" ,
        )

    def set_global_attributes(self, date, monthly=False):
        duration_int = 30 if monthly else 1
        duration_str = 'P1M' if monthly else 'P1D'
        title = "Arctic Sea Ice Age Climate Data Record v2.1",
        if monthly: 
            title = f"{title} (monthly)"
        global_attributes = dict(
            title = title,
            summary = "This climate data record of sea ice age is obtained from coarse resolution ice drift and conentration OSI SAF products. The processing chain features: 1) Lagrangian advection of ice age fractions, 2) Weighted averaging of fractions.",
            iso_topic_category = "oceans,climatology,meteorology,atmosphere",
            keywords = "GCMDSK:Earth Science > Cryosphere > Sea Ice > Sea Ice Concentration, GCMDSK:Earth Science > Oceans > Sea Ice > Sea Ice Concentration, GCMDSK:Earth Science > Climate Indicators > Cryospheric Indicators > Sea Ice Concentration, GCMDSK:Earth Science > Cryosphere > Sea Ice > Sea Ice Motion, GCMDSK:Earth Science > Oceans > Sea Ice > Sea Ice Motion, GCMDSK:Earth Science > Climate Indicators > Cryospheric Indicators > Sea Ice Motion, GCMDLOC:Geographic Region > Northern Hemisphere, GCMDLOC:Vertical Location > Sea Surface, GCMDPROV: CONSORTIA/INSTITUTIONS > NERSC > Nansen Environmental and Remote Sensing Centre",
            keywords_vocabulary = "GCMDSK:GCMD Science Keywords:https://gcmd.earthdata.nasa.gov/kms/concepts/concept_scheme/sciencekeywords, GCMDPROV:GCMD Providers:https://gcmd.earthdata.nasa.gov/kms/concepts/concept_scheme/providers, GCMDLOC:GCMD Locations:https://gcmd.earthdata.nasa.gov/kms/concepts/concept_scheme/locations",
            geospatial_lat_max = 90.00,
            geospatial_lat_min = 17.61,
            geospatial_lon_max = 180.00,
            geospatial_lon_min = -180.00,
            geospatial_vertical_min = 0.00,
            geospatial_vertical_max = 0.00,
            instrument = "SSM/I,SSMIS,AMSR-E,AMSR2",
            platform = "DMSP-F08,DMSP-F10,DMSP-F11,DMSP-F13,DMSP-F14,DMSP-F15,DMSP-F16,DMSP-F17,DMSP-F18,Aqua,GCOM-W1",
            source = self.global_attributes_source,
            time_coverage_start = date.strftime("%Y-%m-%dT00:00:00Z"),
            time_coverage_end = (date + timedelta(duration_int)).strftime("%Y-%m-%dT00:00:00Z"),
            time_coverage_duration = duration_str,
            time_coverage_resolution = duration_str,
            project = "Thickness of Arctic sea ice Reconstructed by Data assimilation and artificial Intelligence Seamlessly (TARDIS)",
            institution = "Nansen Environmental and Remote Sensing Centre",
            creator_name = "Anton Korosov",
            creator_type = "person",
            creator_url = "https://nersc.no/en/ansatt/anton-korosov/",
            creator_email = "anton.korosov@nersc.no",
            creator_institution = "Nansen Environmental and Remote Sensing Centre (NERSC)",
            license = "All intellectual property rights of the Sea Ice Age product belong to NERSC. The use of these products is granted to every user, free of charge. If users wish to use these products, NERSC\'s copyright credit must be shown by displaying the words \'Copyright NERSC\' under each of the products shown. NERSC offers no warranty and accepts no liability in respect of the Sea Ice Age products. NERSC neither commits to nor guarantees the continuity, availability, or quality or suitability for any purpose of, the Sea Ice Age product.",
            references = "https://doi.org/10.5194/tc-12-2073-2018 (Scientific publication)",
            date_created = "2025-07-01T00:00:00Z",
            cdm_data_type = "Grid",
            spatial_resolution = "25.0 km grid spacing",
            algorithm = "lagrangian_sea_ice_age_v2p1",
            geospatial_bounds_crs = "EPSG:6931",
            contributor_name = "Leo Edel, Laurent Bertino",
            contributor_role = "assistant, project leader",
            naming_authority = "NERSC",
            Conventions = "CF-1.7, ACDD-1.3",
            standard_name_vocabulary = "CF Standard Name Table (Version 78, 21 September 2021)",
            product_name = "nersc_arctic_sea_ice_age_climate_data_record",
            product_id = "arctic25km_sea_ice_age_v2p1",
            product_version = "2.1",
            product_status = "released",
            dataset_doi = "10.5281/zenodo.15773501",
            id = "10.5281/zenodo.15773501",
            history = f"Created on 2025-07-01T00:00:00Z by NERSC Sea Ice Age processing script",
        )
        for key, value in global_attributes.items():
            self.setncattr(key, value)

    def set_variable(self, vname, data, dims, atts, dtype=np.float32):
        """
        set variable data and attributes
        Parameters:
        -----------
        vname : str
            name of new variable
        data : numpy.ndarray
            data to set in variable
        dims : list(str)
            list of dimension names for the variable
        atts : dict
            netcdf attributes to set
        dtype : type
            netcdf data type for new variable (eg np.float32 or np.double)
        """
        ncatts = {k:v for k,v in atts.items() if k != '_FillValue'}
        kw = dict(zlib=True)# use compression
        if '_FillValue' in atts:
            # needs to be a keyword for createVariable and of right data type
            kw['fill_value'] = dtype(atts['_FillValue'])
        if 'missing_value' in atts:
            # needs to be of right data type
            ncatts['missing_value'] = dtype(atts['missing_value'])
        dst_var = self.createVariable(vname, dtype, dims, **kw)
        ncatts['grid_mapping'] = self.grid_mapping_variable
        dst_var.setncatts(ncatts)
        dst_var[0] = data

class ExportNetcdf:
    iCDR_start_date = datetime(2021,1,1)
    source_CDR = "Global Sea Ice Drift Climate Data Record Version 1 from the EUMETSAT OSI SAF, \n Sea Ice Concentration Climate Data Record Version 3 from the EUMETSAT OSI SAF"
    source_iCDR = "Daily Low Resolution Sea Ice Displacement from OSI SAF EUMETSAT (OSI-405), \n Sea Ice Concentration Interim Climate Data Record Version 3 from the EUMETSAT OSI SAF"

    def __init__(self, age_grd_dir, dst_root_dir, osisaf_sic_dir, monthly=False, force=False):
        self.age_grd_dir = age_grd_dir
        self.dst_root_dir = dst_root_dir
        self.osisaf_sic_dir = osisaf_sic_dir
        self.monthly = monthly
        self.force = force

    def set_attributes(self):
        self.time_atts = {}
        template_file = f'{self.osisaf_sic_dir}/1991/01/ice_conc_nh_ease2-250_cdr-v3p0_199101011200.nc'
        with Dataset(template_file) as template_ds:
            self.xc = template_ds['xc'][:] * 1000
            self.yc = template_ds['yc'][:] * 1000
            #lon = template_ds['lon'][:]
            #lat = template_ds['lat'][:]
            time_var = template_ds['time']
            for key in time_var.ncattrs():
                self.time_atts[key] = time_var.getncattr(key)
            status_flag = template_ds['status_flag'][0]

        self.status_flag = (status_flag == 1).astype(int)

        self.age_atts = dict(
            standard_name = 'age_of_sea_ice',
            long_name = 'Weighted Average of Sea Ice Age',
            name = 'sea_ice_age',
            comment = 'The weighted average is computed over all available fractions.',
            units = "years",
            valid_min = np.float32(0),
            valid_max = np.float32(20),
            ancillary_variables = "sea_ice_age_uncertainty status_flag",
        )

        self.age_unc_atts = dict(
            long_name = 'Uncertainty in Sea Ice Age',
            name = 'sea_ice_age_uncertainty',
            ancillary_variables = 'status_flag',
            comment = 'Uncertainty is computed based on observational and model errors.',
            units = "years",
            valid_min = np.float32(0),
            valid_max = np.float32(20),
        )

        self.conc_atts = dict(
            long_name = "Concentration of $Numeral$ Year Sea Ice",
            name = "conc_$YEAR$yi",
            units = "1",
            valid_min = np.float32(0),
            valid_max = np.float32(1),
            ancillary_variables = "conc_$YEAR$yi_uncertainty status_flag",
        )

        self.conc_unc_atts = dict(
            long_name = "Uncertainty in concentration of $Numeral$ Year Sea Ice",
            name = "conc_$YEAR$yi_uncertainty",
            units = "1",
            valid_min = np.float32(0),
            valid_max = np.float32(1),
        )

        self.status_atts = dict(
            long_name = "status flag array for sea ice age",
            valid_min = np.byte(1),
            valid_max = np.byte(3),
            grid_mapping = "Lambert_Azimuthal_Grid",
            flag_masks = (np.byte(1), np.byte(2), np.byte(3)),
            flag_meanings = "nominal land invalid",
            flag_descriptions = [
                "flag = 1: Nominal retrieval by the SIA algorithm",
                "flag = 2: Position is over land",
                "flag = 3: Pixel is invalid"],
        )

        self.numerals = [None, 'first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh', 'eighth']

    def __call__(self, age_grd_file):
        age_grd_date = datetime.strptime(os.path.basename(age_grd_file).split('_')[1], '%Y%m%d')
        dst_dir = age_grd_date.strftime(f'{self.dst_root_dir}/%Y')
        os.makedirs(dst_dir, exist_ok=True)
        dst_file = age_grd_date.strftime(f'{dst_dir}/arctic25km_sea_ice_age_v2p1_%Y%m%d.nc')
        if os.path.exists(dst_file) and not self.force:
            return

        time_value = (datetime(age_grd_date.year, age_grd_date.month, age_grd_date.day, 0) - datetime(1978,1,1)).total_seconds()
        time_data = np.array([time_value + 43200], float)
        time_bnds_data = np.array([time_value, time_value+86400], float)[None]
        if self.monthly:
            time_data = np.array([time_value], float)
            time_bnds_data = np.array([time_value, time_value+86400*30], float)[None]

        age_grd = dict(np.load(age_grd_file))
        age_grd_vars = list(age_grd.keys())
        status_flag = self.status_flag.copy()
        status_flag[(self.status_flag == 0) * np.isnan(age_grd['age'])] = 2

        with SeaIceAgeDataset(dst_file, 'w') as ds:
            ds.global_attributes_source = self.source_CDR
            if age_grd_date >= self.iCDR_start_date:
                ds.global_attributes_source = self.source_iCDR
            ds.set_projection_variable()
            ds.set_xy_dims(self.xc, self.yc)
            ds.set_global_attributes(age_grd_date, monthly=self.monthly)
            ds.set_time_variable(time_data, self.time_atts)
            # set_time_bnds_variable
            ds.createDimension('nv', 2)
            tbvar = ds.createVariable('time_bnds', 'f8', ('time', 'nv'), zlib=True)
            tbvar[:] = time_bnds_data

            ds.set_variable('sea_ice_age', age_grd['age'], ('time', 'y', 'x'), self.age_atts, dtype=np.float32)
            ds.set_variable('sea_ice_age_uncertainty', age_grd['age_unc'], ('time', 'y', 'x'), self.age_unc_atts, dtype=np.float32)

            for age_grd_var in sorted(age_grd_vars):
                if 'fraction' in age_grd_var:
                    year = age_grd_var.split('_')[1][0]
                    numeral = self.numerals[int(year)]
                    frac_atts = {}
                    if 'unc' not in age_grd_var:
                        for key, value in self.conc_atts.items():
                            if isinstance(value, str):
                                frac_atts[key] = value.replace('$Numeral$', numeral.title()).replace('$YEAR$', year).replace('$numeral$', numeral)
                            else:
                                frac_atts[key] = value
                    else:
                        for key, value in self.conc_unc_atts.items():
                            if isinstance(value, str):
                                frac_atts[key] = value.replace('$Numeral$', numeral.title()).replace('$YEAR$', year).replace('$numeral$', numeral)
                            else:
                                frac_atts[key] = value
                    ds.set_variable(frac_atts['name'], age_grd[age_grd_var][None]/100., ('time', 'y', 'x'), frac_atts, dtype=np.float32)

            ds.set_variable('status_flag', status_flag + 1, ('time', 'y', 'x'), self.status_atts, dtype=np.byte)
