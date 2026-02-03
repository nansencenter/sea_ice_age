[![DOI](https://zenodo.org/badge/77685919.svg)](https://doi.org/10.5281/zenodo.16743289)

# Sea Ice Age Computation, v2.2

This repository contains code for computing the age of sea ice using satellite data and numerical algorithms. The project is designed to process input datasets, apply age-tracking logic, and output maps or statistics of sea ice age.

## Features

- Ingests satellite-derived sea ice concentration and motion data
- Advects triangular mesh for reducing diffusion
- Outputs age maps in standard formats (e.g., NetCDF, GeoTIFF)
- Modular and extensible codebase

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/yourusername/sea_ice_age.git
cd sea_ice_age
pip install -r requirements.txt
```

## Usage

Create directories that contain OSISAF sea ice concentration (SIC) and sea ice drift (SID)

Run the notebooks in the following order:

01preprocess_sid.ipynb - collocation of SID and SIC products, reprojection on EASE2 grid

02move_mesh.ipynb - advect triangular mesh using SID

03interpolate_sic.ipynb - interpolate SIC onto advected mesh elements

04find_autumn_minima.ipynb - find MYI concentration for 15 September each year

05compute_age.ipynb - propagate MYI concentration using the advected mesh and use the propagated MYI concentration from different start years to compute ice age fractions

07uncert_observations.ipynb - compute uncertainty of sea ice drift and sea ice concentration

08uncert_integrate.ipynb - integrate drift and concentration uncertainties for the advected MYI concentration

09uncert_age.ipynb - compute total uncertainty

11grid_age.ipynb - interpolate the age fractions from the advected triangular mesh onto a grid

12netcdf_age.ipynb - create netCDFs

12netcdf_age_cmems.ipynb - create netCDFs for CMEMS

Directories `iabp_validation` and `figures` contain notebooks for the ESSD manuscript. The directory `nextsim` contains code to compute the ice age from neXtSIM mesh files.

## Data Sources

- OSISAF Sea Ice Concentration CDR
- OSISAF Sea Ice Drift CDR and iCDR

## Output

- Gridded maps of sea ice age fractions, and weighted average sea ice age

## License

See [LICENSE.txt](LICENSE) for details.

## Contact

For questions or contributions, open an issue or contact the maintainer.

## Citation

Korosov, A., Edel, L., Regan, H., Lavergne, T., Aaboe, S., and Down, E. J.: A climate data record of sea ice age using Lagrangian advection of a triangular mesh, Earth Syst. Sci. Data, 18, 721–740, https://doi.org/10.5194/essd-18-721-2026, 2026.

Korosov, A. A., Rampal, P., Pedersen, L. T., Saldo, R., Ye, Y., Heygster, G., Lavergne, T., Aaboe, S., and Girard-Ardhuin, F.: A new tracking algorithm for sea ice age distribution estimation, The Cryosphere, 12, 2073–2085, https://doi.org/10.5194/tc-12-2073-2018, 2018.

