# Original sea-ice-age output contract

The CMIP6 workflow preserves the final variable names and meanings used by the
observational `ExportNetcdf` implementation wherever the model inputs provide
the required information.

## Always exported

- `sea_ice_age`: weighted average sea-ice age in years.
- `conc_1yi` through `conc_6yi`: concentrations of the available age classes,
  stored as fractions from 0 to 1.
- `status_flag`: `1` nominal, `2` land, `3` invalid.
- `time` and `time_bnds`.
- `x` and `y` projection coordinates.
- `Lambert_Azimuthal_Equal_Area` grid-mapping variable.
- `longitude` and `latitude` when enabled in configuration.

Age classes that do not yet exist at the start of a run are omitted rather
than filled with invented data. A multi-year run progressively produces the
same set of age-class variables as the original algorithm, up to the configured
six-class limit used by `05compute_age.ipynb`.

## Conditional uncertainty variables

The original observational product can additionally contain
`sea_ice_age_uncertainty` and `conc_Nyi_uncertainty`. Those variables are
calculated by notebooks 07, 08, and 09 from observational SIC and drift
uncertainties. The current CMCC input has no corresponding uncertainty fields.
They are therefore disabled rather than populated with artificial zeros or
missing values. Their names remain reserved and the exporter can add them when
scientifically defensible uncertainty inputs are configured.

## Model-specific metadata

Satellite platform, instrument, source, DOI, and observational uncertainty
statements are replaced by CMCC/CMIP6 provenance. The numerical variable
definitions, units, status meanings, projection, and age-class naming remain
compatible with the original output.
