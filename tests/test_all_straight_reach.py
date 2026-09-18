import unittest

import numpy as np
import pandas as pd
from shapely.geometry import LineString

from pipeline.inflection_points import inflection_points_curve
from pipeline.get_orthogonals import build_orthogonal_lines


class AllStraightReachTests(unittest.TestCase):
    def test_consecutive_straight_candidates_retain_whole_reach(self):
        line = LineString([(0, 0), (100, 2), (200, -2), (300, 1), (400, 0)])
        nodes = pd.DataFrame({
            'reach_id': [1] * 5, 'linePos': np.linspace(0, line.length, 5),
            'width': [5.] * 5, 'max_width': [5.] * 5,
            'dist_out': [500., 400., 300., 200., 100.], 'facc': [10.] * 5,
        })
        reach = pd.DataFrame({'reach_id': [1], 'width': [5.],
                              'max_width': [5.], 'combined_reach_id': [1]})
        result = inflection_points_curve(line, reach, nodes)
        self.assertEqual(len(result[8]), 1)
        np.testing.assert_array_equal(result[8][0].coords, line.coords)
        np.testing.assert_array_equal(result[9][0].coords, [line.coords[0], line.coords[-1]])
        self.assertEqual(result[0], line.length / 400.)
        self.assertEqual(result[13][0], line.length)
        orthogonals = build_orthogonal_lines(line, result[5], result[6], result[4],
                                            result[9], result[10], 50)
        for geometry in (*orthogonals.line_out, *orthogonals.line_inn):
            self.assertTrue(geometry.is_valid)
            self.assertGreater(geometry.length, 0)


if __name__ == '__main__':
    unittest.main()
