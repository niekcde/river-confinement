import unittest

import numpy as np
from shapely.geometry import LineString

from pipeline.line_functions import get_points_along_linestring


def scalar_reference(line, spacing=20):
    distances = np.linspace(0, line.length, int(line.length / spacing) + 1)
    return [list(line.interpolate(distance).coords[0]) for distance in distances]


class LineSamplingTests(unittest.TestCase):
    def test_exact_scalar_equivalence(self):
        lines = [
            LineString([(0, 0), (100, 100)]),
            LineString([(100, 100), (0, 0)]),
            LineString([(0, 0), (20, 80), (100, 30)]),
            LineString([(0, 0), (100, 100), (0, 100), (100, 0)]),
            LineString([(0, 0), (100, 0), (0, 0)]),
            LineString([(5, 5), (5, 5)]),
            LineString([(0, 0, 1), (20, 80, 4), (100, 30, 2)]),
        ]
        for line in lines:
            for spacing in (1, 10, 20, 1000):
                with self.subTest(line=line.wkt, spacing=spacing):
                    expected = scalar_reference(line, spacing)
                    actual = get_points_along_linestring(line, spacing)
                    self.assertIsInstance(actual, list)
                    self.assertTrue(all(isinstance(point, list) for point in actual))
                    np.testing.assert_array_equal(actual, expected, strict=True)


if __name__ == "__main__":
    unittest.main()
