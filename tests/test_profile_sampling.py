import unittest

import numpy as np
import xarray as xr
from shapely.geometry import LineString

from pipeline.extract_slope_along_raster_line import extract_slope_along_raster_line


def scalar_reference(raster, line, samples=400):
    """Original sampler, retained to verify numerical compatibility."""
    positions = np.arange(samples, dtype=float) / samples - 1.0
    points = [line.interpolate(p, normalized=True) for p in positions]
    xs = xr.DataArray([p.x for p in points], dims="points")
    ys = xr.DataArray([p.y for p in points], dims="points")
    distances = [line.project(p) for p in points]
    profile = raster.sel(x=xs, y=ys, method="nearest").data
    if len(profile) == 1:
        profile = profile[0]
    return distances, profile


class ProfileSamplingTests(unittest.TestCase):
    def test_exact_scalar_equivalence(self):
        raster = xr.DataArray(
            np.arange(121, dtype=float).reshape(11, 11),
            dims=("y", "x"), coords={"y": np.arange(10, -1, -1), "x": np.arange(11)},
        )
        raster.values[5, 5] = np.nan
        lines = [
            LineString([(0, 0), (10, 10)]),
            LineString([(10, 10), (0, 0)]),
            LineString([(0, 0), (2, 8), (10, 3)]),
            LineString([(0, 0), (10, 10), (0, 10), (10, 0)]),
            LineString([(0, 0), (10, 0), (0, 0)]),
            LineString([(5, 5), (5, 5)]),
        ]
        for source in (raster, raster.expand_dims(band=[1])):
            for line in lines:
                for samples in (1, 2, 17, 400):
                    with self.subTest(banded=source.ndim == 3, line=line.wkt, samples=samples):
                        expected = scalar_reference(source, line, samples)
                        actual = extract_slope_along_raster_line(source, line, samples)
                        self.assertIsInstance(actual[0], list)
                        np.testing.assert_array_equal(actual[0], expected[0], strict=True)
                        np.testing.assert_array_equal(actual[1], expected[1], strict=True)


if __name__ == "__main__":
    unittest.main()
