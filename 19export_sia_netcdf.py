import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call

from nansat import *

from iceagelib import *
from cmocean import cm

### EXPORT SIA, MYI, SIF into netcdf
sif_names_selector = {
0: 'first',
1: 'second',
2: 'third',
3: 'fourth',
4: 'fifth',
5: 'sixth',
6: 'seventh',
7: 'eighth',
8: 'ninth',
9: 'tenth',
}

call('rm test.nc', shell=True)

# source domain
osi_nsr = NSR('+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45')
osi_sic_dom = Domain(osi_nsr, '-te -3850000 -5350000 3750000 5850000 -tr 20000 20000')

osi_sia_dir = '/files/sea_ice_age/nersc_osi_fv6_2017_conc/'

out_dir = osi_sia_dir + 'netdcf/'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

for yy in [2012, 2013, 2014, 2015, 2016, 2017]:
    # export to netcdf
    idir = osi_sia_dir + 'sia/'
    ifiles = sorted(glob.glob(idir + '*%d-??-??_sia.npz' % yy))
    for ifile in ifiles:
        print ifile
        idate = parse(os.path.basename(ifile).split('_')[0])
        ofile = os.path.join(out_dir, 'arctic20km_sia_%s.nc' % idate.strftime('%Y%m%d'))
        if os.path.exists(ofile):
            continue

        time_coverage_start = idate.isoformat()
        time_coverage_end = (idate+dt.timedelta(1)).isoformat()
        
        sif = np.load(ifile)['sif'][:, ::2, ::2]
        sia = np.load(ifile)['sia'][::2, ::2]
        myi = np.load(ifile)['myi'][::2, ::2]
        wat = sif.sum(axis=0) == 0
        lnd = np.isnan(sia)
        sia[wat] = np.nan
        sif[:, lnd] = np.nan
        myi[lnd] = np.nan
        # bitwise flag
        # 0bit - valid product
        # 1bit - open water
        # 2bit - land
        flg = np.ones(lnd.shape, np.uint8)
        flg[wat] = 2
        flg[lnd] = 4
        
        bands = [sia, myi, flg]
        parameters=[{'name':'sia'}, {'name':'conc_myi'}, {'name':'status_flag'}]
        bands_dict = {
        'sia': {'type'     : '>f4',
                'standard_name' : 'age_of_sea_ice',
                'long_name': 'Weighted Average of Sea Ice Age',
                'ancillary_variables': 'status_flag',
                'units' : 'year'},
        'conc_myi': {'type'     : '>f4',
                'standard_name' : 'multi_year_sea_ice_area_fraction',
                'long_name': 'Concentration of Multiyear Sea Ice',
                'units' : '1'},
        'status_flag': {'type'     : '>i1',
                'standard_name' : 'age_of_sea_ice status_fag',
                'long_name': 'Status Flag of Sea Ice Age',
                'flag_values': [1, 2, 4],
                'flag_meanings': 'nominal_quality over_open_water over_land'},
        }
        metadata = {
            'time_coverage_start': time_coverage_start,
            'time_coverage_end': time_coverage_end,
            'license' : 'CC BY 4.0',
            'comment' : 'Supplementary material for (Korosov et al., A new Tracking Algrotihm for Sea Ice Age Distribution Estimation, The Cryosphere, submitted)',
            'Metadata_Conventions' : 'Unidata Dataset Discovery v1.0',
            'keywords' : 'EARTH SCIENCE > CRYOSPHERE > SEA ICE > SEA ICE AGE',
            'keywords_vocabulary' : 'GCMD Science Keywords',
            'id' : 'arctic20km_sia_%s' % idate.strftime('%Y%m%d'),
            'time_coverage_resolution' : 'P1D',
            'acknowledgment' : 'European SPace Agency, Sea Ice Climate Change Initiative, Sea Ice Age Option',
            'title' : 'Arctic Sea Ice Age (SIA)',
            'naming_authority' : 'nersc.no',
            'processing_level' : '4',
            'date_issued' : '2017-11-01',
            'institution' : 'Nansen Environmental and Remote Sensing Centre',
            'summary' : 'Sea ice age in the Arctic derived using the new tracking algorithm for sea ice age distribution derived from ice motion and concentration',
            'cdm_data_type' : 'Grid',
            'project' : 'ESA SICCI SIA',
            'publisher_name' : 'NERSC (Anton Korosov)',
            'publisher_email' : 'anton.korosov@nersc.no',
            'publisher_url' : 'http//www.nersc.no',
        }
        
        for sif_index, sif_conc in enumerate(sif):
            if sif_index == 0:
                sif_name = 'conc_fyi'
            else:
                sif_name = 'conc_%dyi' % (sif_index + 1)
            sif_std_name = '%s_year_sea_ice_area_fraction' % sif_names_selector[sif_index]
            sif_lng_name = 'Concentration of %s Year Sea Ice' % sif_names_selector[sif_index].capitalize()
            
            parameters.append({'name': sif_name})
            bands.append(sif_conc)
            bands_dict[sif_name] = {
                    'type'     : '>f4',
                    'standard_name' : sif_std_name,
                    'long_name' : sif_lng_name,
                    'units' : '1'}

        n = Nansat(domain=osi_sic_dom)
        n.add_bands(bands, parameters)
        n.set_metadata('time_coverage_start', time_coverage_start)
        n.set_metadata('time_coverage_end', time_coverage_end)
        n.export2thredds(ofile, bands_dict, metadata)
        n = None
        del n
        #call('ncdump -h test.nc', shell=True)
