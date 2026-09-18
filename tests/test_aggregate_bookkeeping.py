import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import xarray as xr

from pipeline.get_orthogonals import _mean_valid_profile_height
from pipeline.run_confinement_values import concat_nc_conf_files, _prepare_step5_netcdf_dataframe


class AggregateBookkeepingTests(unittest.TestCase):
    def test_global_indices_names_and_missing_elevation(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            names = ['A', 'B', 'Dintel; Hellegat; Schelde-Rijnkanaal', 'Rivière très longue']
            for batch in range(2):
                ds = xr.Dataset({
                    'river_name': ('index', names[batch * 2:batch * 2 + 2]),
                    'bendHeight': ('index', [10., -9999.] if batch == 0 else [99999., -5.]),
                    'slope_out': ('index', [-0.01, 0.03]),
                }, coords={'index': [0, 1]})
                ds.to_netcdf(folder / f'eu_{batch:02}_50_02_conf.nc',
                             encoding={'river_name': {'dtype': 'S1'}})
            with patch('pipeline.run_confinement_values._resolve_step6_paths',
                       return_value={'single_values_dir': folder}):
                result = concat_nc_conf_files()
            with xr.open_dataset(result) as ds:
                self.assertEqual(ds.river_name.values.tolist(), names)
                np.testing.assert_array_equal(ds['index'], [0, 1, 2, 3])
                np.testing.assert_array_equal(ds.source_index, [0, 1, 0, 1])
                np.testing.assert_array_equal(ds.bendHeight, [10., np.nan, np.nan, -5.])
                np.testing.assert_array_equal(ds.slope_out, [-0.01, 0.03, -0.01, 0.03])

    def test_step5_normalizes_sentinels_not_real_negative_heights(self):
        df = pd.DataFrame({'bendHeight': [-9999., 99999., -5., 0., 20.]})
        result = _prepare_step5_netcdf_dataframe(df)
        np.testing.assert_array_equal(result.bendHeight, [np.nan, np.nan, -5., 0., 20.])

    def test_height_mean_ignores_missing_preserves_negative_elevation(self):
        self.assertEqual(_mean_valid_profile_height([-9999, 99999, np.nan, np.inf, -10, 20], -9999), 5.)
        self.assertTrue(np.isnan(_mean_valid_profile_height([-9999, 99999, np.nan], -9999)))


if __name__ == '__main__':
    unittest.main()
