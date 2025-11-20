from datetime import datetime, timedelta
import os

from netCDF4 import Dataset
from geodataset.geodataset import GeoDatasetWrite
import numpy as np
import pyproj

from lmsiage.utils import IrregularGridInterpolator
from lmsiage.mesh_file import MeshFile
from lmsiage.zarr_index import get_file_arrays

class GridAge:
    def __init__(self, mesh_dir, age_dir, unc_dir, grid_dir, xgrd, ygrd, mask, force=False):
        self.mesh_dir = mesh_dir
        self.age_dir = age_dir
        self.unc_dir = unc_dir
        self.grid_dir = grid_dir
        self.xgrd = xgrd
        self.ygrd = ygrd
        self.mask = mask
        self.force = force

    def load_data(self, mesh_file, age_file, unc_file):
        """ Load mesh, age, and uncertainty data from files."""
        self.x, self.y, self.t = MeshFile(mesh_file).load(['x', 'y', 't'], as_dict=False)
        self.age, self.fracs, = MeshFile(age_file).load(['age', 'f'], as_dict=False)
        self.unc_age, self.unc_fracs, = MeshFile(unc_file).load(['unc_age', 'unc_fracs'], as_dict=False)

    def __call__(self, unc_file):
        """ Load, interpolate, save grid """
        date = datetime.strptime(os.path.basename(unc_file), 'unc_%Y%m%d.zip')
        mesh_file = date.strftime(f'{self.mesh_dir}/%Y/mesh_%Y%m%d.zip')
        age_file = date.strftime(f'{self.age_dir}/%Y/age_%Y%m%d.zip')
        grid_dir = date.strftime(f'{self.grid_dir}/%Y')
        grid_file = date.strftime(f'{grid_dir}/grid_%Y%m%d.zip')
        
        file_arrays = get_file_arrays(grid_file)
        if 'age' in file_arrays and not self.force:
            return

        self.load_data(mesh_file, age_file, unc_file)
        igi = IrregularGridInterpolator(self.x, self.y, self.xgrd, self.ygrd, self.t)
        src_data = np.vstack([self.age[None], self.fracs, self.unc_age[None], self.unc_fracs])
        frac_num = len(self.fracs)
        save_data_names = (
            ['age'] + 
            [f'sic_{frac_num - i}yi' for i in range(frac_num)] +
            ['unc_age'] + 
            [f'unc_{frac_num - i}yi' for i in range(frac_num)]
        )
        dst_data = {}
        for d, name in zip(src_data, save_data_names):
            dgrd = igi.interp_field(d)
            dgrd[self.mask == 0] = np.nan
            dst_data[name] = dgrd.astype(np.float16)

        os.makedirs(grid_dir, exist_ok=True)
        mf = MeshFile(grid_file)
        mf.save(dst_data)


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
            acknowledgment = " Research Council of Norway (project 'TARDIS', no. 325241) and the European Space Agency (project 'CCI SAGE', no. 4000147560/25/I-LR",
            comment = 'No comments',
            geospatial_bounds = "POLYGON((-5387500 -5387500, 5387500 -5387500, 5387500 5387500, -5387500 5387500, -5387500 -5387500))",
            processing_level = "Level 4",
            publisher_name = "Anton Korosov",
            publisher_url = "https://nersc.no/en/ansatt/anton-korosov/",
            publisher_email = "anton.korosov@nersc.no",
        )
        for key, value in global_attributes.items():
            self.setncattr(key, value)

    def set_variable(self, vname, data, dims, atts, dtype=np.float32, pack=None):
        """
        pack: optional dict with keys:
              dtype (np.int16/np.uint8/etc), scale_factor (float), add_offset (float, default 0),
              _FillValue (int), missing_value (int, optional)
        """
        if pack is not None:
            # add packing attrs to metadata (valid_* stay in real units)
            scale_factor = np.float32(pack.get('scale_factor', 1.0))
            add_offset = np.float32(pack.get('add_offset', 0.0))

            atts = {**atts,
                    'scale_factor': scale_factor,
                    'add_offset': add_offset,
                    '_FillValue': np.iinfo(dtype).min}

            atts['valid_max'] = dtype((atts['valid_max'] - add_offset)/ scale_factor)
            atts['valid_min'] = dtype((atts['valid_min'] - add_offset)/ scale_factor)

            fill_value_scaled = np.iinfo(dtype).min * scale_factor + add_offset
            data[~np.isfinite(data)] = fill_value_scaled

        ncatts = {k: v for k, v in atts.items() if k != '_FillValue'}
        kw = dict(zlib=True)
        if '_FillValue' in atts:
            kw['fill_value'] = dtype(atts['_FillValue'])
        if 'missing_value' in atts:
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
            coverage_content_type = "modelResult",
        )

        self.age_unc_atts = dict(
            long_name = 'Uncertainty in Sea Ice Age',
            standard_name = 'age_of_sea_ice standard_error',
            name = 'sea_ice_age_uncertainty',
            ancillary_variables = 'status_flag',
            comment = 'Uncertainty is computed based on observational and model errors.',
            units = "years",
            valid_min = np.float32(0),
            valid_max = np.float32(21),
            coverage_content_type = "modelResult",
        )

        self.conc_atts = dict(
            long_name = "Concentration of $Numeral$ Year Sea Ice",
            name = "conc_$YEAR$yi",
            units = "1",
            valid_min = np.float32(0),
            valid_max = np.float32(1.1),
            ancillary_variables = "conc_$YEAR$yi_uncertainty status_flag",
            coverage_content_type = "modelResult",
            standard_name='sea_ice_area_fraction',
        )

        self.conc_unc_atts = dict(
            long_name = "Uncertainty in concentration of $Numeral$ Year Sea Ice",
            standard_name = 'sea_ice_area_fraction standard_error',            
            name = "conc_$YEAR$yi_uncertainty",
            units = "1",
            valid_min = np.float32(0),
            valid_max = np.float32(1.1),
            coverage_content_type = "modelResult",
        )

        self.status_atts = dict(
            standard_name = 'age_of_sea_ice status_flag',            
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
            coverage_content_type = "qualityInformation",
        )

        self.numerals = [None, 'first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh', 'eighth']

    def __call__(self, age_grd_file):
        age_grd_date = datetime.strptime(os.path.basename(age_grd_file), 'grid_%Y%m%d.zip')
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

        age_grd = MeshFile(age_grd_file).load()
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
            unc_age_masked = np.ma.masked_array(age_grd['unc_age'], mask=np.isnan(age_grd['age']))
            unc_age_masked = np.clip(unc_age_masked, 0., 20.)
            ds.set_variable('sea_ice_age_uncertainty', unc_age_masked, ('time', 'y', 'x'), self.age_unc_atts, dtype=np.int16, pack=dict(scale_factor=0.001, add_offset=0.0))

            for age_grd_var in sorted(age_grd_vars):
                if age_grd_var.endswith('yi'):
                    if age_grd_var.startswith('sic_'):
                        use_attributes = self.conc_atts.items()
                    elif age_grd_var.startswith('unc_'):
                        use_attributes = self.conc_unc_atts.items()

                    dtype = np.int16
                    pack = dict(scale_factor=0.001, add_offset=0.0)
                    output_array = age_grd[age_grd_var][None]/100.
                    output_array = np.clip(output_array, 0., 1.)
                    output_array = np.ma.masked_array(output_array, mask=np.isnan(age_grd[age_grd_var]))

                    year = age_grd_var.split('_')[1][0]
                    numeral = self.numerals[int(year)]
                    frac_atts = {}
                    for key, value in use_attributes:
                        if isinstance(value, str):
                            frac_atts[key] = value.replace('$Numeral$', numeral.title()).replace('$YEAR$', year).replace('$numeral$', numeral)
                        else:
                            frac_atts[key] = value
                    ds.set_variable(frac_atts['name'], output_array, ('time', 'y', 'x'), frac_atts, dtype=dtype, pack=pack)

            ds.set_variable('status_flag', status_flag + 1, ('time', 'y', 'x'), self.status_atts, dtype=np.byte)
