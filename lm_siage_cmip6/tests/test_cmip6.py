from pathlib import Path

import pytest

from lm_siage_cmip6.cmip6 import parse_cmip6_filename
from lm_siage_cmip6.pipeline import annual_age_filename


def test_parse_cmip6_filename_with_hyphenated_model():
    item = parse_cmip6_filename(
        "siconc_SIday_NorESM2-MM_historical_r1i1p1f1_gn_18500101-18591231.nc"
    )
    assert item.variable == "siconc"
    assert item.model == "NorESM2-MM"
    assert item.experiment == "historical"
    assert item.variant == "r1i1p1f1"
    assert item.start.year == 1850
    assert item.stop.year == 1859


def test_reject_non_cmip6_output_name():
    with pytest.raises(ValueError):
        parse_cmip6_filename(Path("RemappedInput_NorESM2-MM_1850_LAEA.nc"))


def test_annual_age_filename():
    config = {
        "dataset": {
            "model": "NorESM2-MM",
            "experiment": "historical",
            "variant": "r1i1p1f1",
        },
        "output": {
            "filename_template": "SeaIceAge_{model}_{experiment}_{variant}_{year}.nc"
        },
    }
    assert annual_age_filename(config, 1853) == (
        "SeaIceAge_NorESM2-MM_historical_r1i1p1f1_1853.nc"
    )
