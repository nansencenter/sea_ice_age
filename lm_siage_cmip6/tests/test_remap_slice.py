from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


SCRIPT = Path(__file__).parents[1] / "scripts" / "00_remap_model_data.py"
SPEC = spec_from_file_location("remap_model_data", SCRIPT)
MODULE = module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_time_slice_is_rebased_to_zero(tmp_path):
    source = tmp_path / "source.nc"
    output = tmp_path / "output.nc"
    with Dataset(source, "w") as dataset:
        dataset.createDimension("time", None)
        dataset.createDimension("x", 2)
        time = dataset.createVariable("time", "i4", ("time",))
        values = dataset.createVariable("values", "f4", ("time", "x"))
        time[:] = np.arange(8)
        values[:] = np.arange(16).reshape(8, 2)

    MODULE.copy_time_slice(source, output, 3, 6)

    with Dataset(output) as dataset:
        assert len(dataset.dimensions["time"]) == 3
        np.testing.assert_array_equal(dataset.variables["time"][:], [3, 4, 5])
        np.testing.assert_array_equal(
            dataset.variables["values"][:], np.arange(16).reshape(8, 2)[3:6]
        )


def test_static_status_is_repeated_across_time(tmp_path):
    output = tmp_path / "remapped.nc"
    status = np.array([[0, 1], [2, 0]], dtype=np.int8)
    with Dataset(output, "w") as dataset:
        dataset.createDimension("time", 3)
        dataset.createDimension("y", 2)
        dataset.createDimension("x", 2)
        variable = dataset.createVariable(
            "status_flag", "i1", ("time", "y", "x"), fill_value=np.int8(-127)
        )
        variable[0] = status

    MODULE.repeat_static_status(output)

    with Dataset(output) as dataset:
        assert int(dataset.static_status_repeated) == 1
        for index in range(3):
            np.testing.assert_array_equal(dataset.variables["status_flag"][index], status)
