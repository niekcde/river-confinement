import unittest

import numpy as np
import pandas as pd
from shapely.geometry import LineString

from pipeline.line_functions import get_bend_width


class BendWidthFallbackTests(unittest.TestCase):
    def setUp(self):
        self.line = LineString([(0, 0), (100, 0)])
        self.nodes = pd.DataFrame({
            'reach_id': np.array([1, 2, 3], dtype=np.int64),
            'linePos': [0., 60., 100.],
            'width': [10., 20., 30.], 'max_width': [15., 25., 35.],
        })
        self.reaches = self.nodes[['reach_id', 'width', 'max_width']].copy()

    def test_nearest_to_absolute_midpoint(self):
        # Midpoint 75: node 60 is nearest; old code incorrectly chose node 0.
        self.assertEqual(get_bend_width(self.line, LineString([(70, 0), (80, 0)]),
                                        self.nodes, self.reaches), (20., 25.))

    def test_unsorted_nodes_and_ties_preserve_row_order(self):
        nodes = self.nodes.iloc[[2, 0, 1]]
        self.assertEqual(get_bend_width(self.line, LineString([(70, 0), (90, 0)]),
                                        nodes, self.reaches), (30., 35.))

    def test_interior_nodes_still_averaged(self):
        self.assertEqual(get_bend_width(self.line, LineString([(0, 0), (100, 0)]),
                                        self.nodes, self.reaches), (20., 25.))

    def test_missing_nearest_width_is_skipped(self):
        self.nodes.loc[1, ['width', 'max_width']] = np.nan
        self.assertEqual(get_bend_width(self.line, LineString([(70, 0), (80, 0)]),
                                        self.nodes, self.reaches), (30., 35.))

    def test_invalid_widths_are_skipped(self):
        for column in ('width', 'max_width'):
            for value in (np.nan, np.inf, 0, -1):
                with self.subTest(column=column, value=value):
                    nodes = self.nodes.copy()
                    nodes.loc[1, column] = value
                    self.assertEqual(get_bend_width(self.line, LineString([(70, 0), (80, 0)]),
                                                    nodes, self.reaches), (30., 35.))

    def test_no_valid_node_has_explicit_error(self):
        self.nodes['width'] = np.nan
        with self.assertRaisesRegex(ValueError, 'No node with finite positive'):
            get_bend_width(self.line, LineString([(70, 0), (80, 0)]), self.nodes, self.reaches)


if __name__ == '__main__':
    unittest.main()
