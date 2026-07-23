from datetime import datetime
from types import SimpleNamespace

from netCDF4 import Dataset, date2num
import numpy as np

from lm_siage_cmip6.catalog import InputCatalog
from lm_siage_cmip6.pipeline import (
    current_experiment_start,
    enforce_cumulative_nesting,
    latest_complete_output,
    propagate_element_field,
    select_recent_cohort_names,
)


def test_identity_propagation():
    source = np.array([20.0, 80.0])
    mapping = np.array([[0, 0], [1, 1]])
    weights = np.ones(2)
    ratio = np.ones(2)
    result = propagate_element_field(source, mapping, weights, ratio)
    np.testing.assert_allclose(result, source)


def test_propagation_is_capped_by_observation():
    result = propagate_element_field(
        np.array([80.0]), np.array([[0, 0]]), np.ones(1), np.ones(1),
        observed_limit=np.array([50.0]),
    )
    np.testing.assert_allclose(result, [50.0])


def test_area_ratio_clip_is_optional():
    source = np.array([50.0])
    mapping = np.array([[0, 0]])
    weights = np.ones(1)
    ratio = np.array([2.0])
    unclipped = propagate_element_field(source, mapping, weights, ratio)
    clipped = propagate_element_field(
        source, mapping, weights, ratio, area_ratio_clip=(0.9, 1.1)
    )
    np.testing.assert_allclose(unclipped, [25.0])
    np.testing.assert_allclose(clipped, [50.0 / 1.1])


def test_cumulative_cohorts_are_nested_oldest_to_newest():
    nested = enforce_cumulative_nesting(
        [np.array([70.0, 20.0]), np.array([50.0, 30.0])],
        np.array([80.0, 25.0]),
    )

    np.testing.assert_allclose(nested[0], [50.0, 20.0])
    np.testing.assert_allclose(nested[1], [50.0, 25.0])
    fractions = np.diff(
        np.asarray([np.zeros(2)] + nested + [np.array([80.0, 25.0])]), axis=0
    )
    assert np.all(fractions >= 0.0)


def test_only_six_newest_annual_cohorts_are_selected():
    data = {
        f"sic{year}0930": np.array([year], dtype=float)
        for year in range(1850, 1857)
    }
    assert select_recent_cohort_names(data, 6) == [
        f"sic{year}0930" for year in range(1851, 1857)
    ]


def _write_complete_output(path, date_key):
    date = datetime.strptime(date_key, "%Y%m%d")
    units = "days since 1850-01-01"
    with Dataset(path, "w") as dataset:
        dataset.createDimension("time", 1)
        time = dataset.createVariable("time", "f8", ("time",))
        time.units = units
        time.calendar = "noleap"
        time[:] = [date2num(date, units, calendar="noleap")]
        dataset.createVariable("sea_ice_age", "f4", ("time",))[:] = [1.0]
        dataset.createVariable("status_flag", "i1", ("time",))[:] = [1]
        dataset.processing_complete = 1


def test_restart_prefers_current_scenario_and_never_uses_another_scenario(tmp_path):
    output = tmp_path / "output" / "NorESM2-MM"
    output.mkdir(parents=True)
    _write_complete_output(
        output / "SeaIceAge_NorESM2-MM_historical_r1i1p1f1_2014.nc",
        "20141231",
    )
    wrong_scenario = output / "SeaIceAge_NorESM2-MM_ssp585_r1i1p1f1_2100.nc"
    _write_complete_output(wrong_scenario, "21001231")
    current_scenario = output / "SeaIceAge_NorESM2-MM_ssp126_r1i1p1f1_2020.nc"
    _write_complete_output(current_scenario, "20201231")
    config = {
        "repo_root": tmp_path,
        "output_dir": tmp_path / "output",
        "dataset": {
            "model": "NorESM2-MM", "experiment": "ssp126", "variant": "r1i1p1f1"
        },
        "restart": {"parent_experiment": "historical"},
    }
    catalog = SimpleNamespace(records=[
        SimpleNamespace(key="20141231"),
        SimpleNamespace(key="20201231"),
        SimpleNamespace(key="21001231"),
    ])

    checkpoint = latest_complete_output(config, catalog)
    assert checkpoint.path == current_scenario
    current_scenario.unlink()
    checkpoint = latest_complete_output(config, catalog)
    assert checkpoint.experiment == "historical"
    assert checkpoint.date_key == "20141231"


def test_current_experiment_start_excludes_bridge(tmp_path):
    bridge = (tmp_path / "bridge.nc").resolve()
    config = {
        "repo_root": tmp_path,
        "restart": {"bridge_input_file": bridge},
    }
    catalog = SimpleNamespace(records=[
        SimpleNamespace(key="20141231", path=bridge),
        SimpleNamespace(key="20150101", path=(tmp_path / "scenario.nc").resolve()),
    ])
    assert current_experiment_start(config, catalog) == "20150101"


def _write_catalog_time(path, dates):
    units = "days since 1850-01-01"
    with Dataset(path, "w") as dataset:
        dataset.createDimension("time", len(dates))
        time = dataset.createVariable("time", "f8", ("time",))
        time.units = units
        time.calendar = "noleap"
        time[:] = date2num(dates, units, calendar="noleap")


def test_catalog_prepends_one_day_historical_bridge(tmp_path):
    bridge = tmp_path / "bridge.nc"
    scenario = tmp_path / "scenario.nc"
    _write_catalog_time(bridge, [datetime(2014, 12, 31)])
    _write_catalog_time(
        scenario, [datetime(2015, 1, 1), datetime(2015, 1, 2)]
    )
    config = {
        "repo_root": tmp_path,
        "input_file": scenario,
        "input_pattern": scenario,
        "restart": {"bridge_input_file": bridge},
    }
    catalog = InputCatalog.from_config(config)
    assert [record.key for record in catalog.records] == [
        "20141231", "20150101", "20150102"
    ]
