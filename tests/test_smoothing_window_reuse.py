import contextlib
import io
import unittest
from unittest.mock import patch

import numpy as np
import scipy as sc
import shapely
from shapely.geometry import LineString

from pipeline import smoothing


def original_filter(line, w):
    if w > len(line.coords):
        w = len(line.coords)
    if w == 0:
        w = 2
    w = int(w)
    x = sc.signal.savgol_filter(line.xy[0], w, 1)
    y = sc.signal.savgol_filter(line.xy[1], w, 1)
    return LineString(list(zip(x, y)))


def original_smoothing(line, w, width, seg=1, simp=0.1,
                       simplify_line=True, id=9999):
    line = line.segmentize(seg)
    if simplify_line:
        line = line.simplify(simp, preserve_topology=True)
    result = original_filter(line, w)
    distance = shapely.hausdorff_distance(line, result)
    while distance > width * 0.5:
        w *= 0.95
        if w < 2:
            print(f'Smoothing broken at window size 2 ({w/0.95}) ({id})')
            break
        result = original_filter(line, w)
        distance = shapely.hausdorff_distance(line, result)
    return result


class SmoothingWindowTests(unittest.TestCase):
    def test_coordinates_and_warnings_match_exactly(self):
        lines = [
            LineString([(0, 0), (30, 0)]),
            LineString([(0, 0), (10, 20), (20, 0), (30, 20)]),
            LineString([(0, 0), (20, 20), (0, 20), (20, 0)]),
        ]
        for line in lines:
            for window in (0, 2.09, 10.8, 200):
                for width in (0, 0.1, 100):
                    for simplify in (True, False):
                        with self.subTest(window=window, width=width, simplify=simplify):
                            old_log, new_log = io.StringIO(), io.StringIO()
                            with contextlib.redirect_stdout(old_log):
                                expected = original_smoothing(line, window, width, simplify_line=simplify)
                            with contextlib.redirect_stdout(new_log):
                                actual = smoothing.SG_smoothing(line, window, width, simplify_line=simplify)
                            np.testing.assert_array_equal(actual.coords, expected.coords, strict=True)
                            self.assertEqual(new_log.getvalue(), old_log.getvalue())

    def test_each_effective_window_is_computed_once(self):
        line = LineString([(0, 0), (10, 20), (20, 0), (30, 20)])
        with patch.object(smoothing, 'apply_smoothing', wraps=smoothing.apply_smoothing) as calls:
            with contextlib.redirect_stdout(io.StringIO()):
                smoothing.SG_smoothing(line, 200, 0, seg=100, simplify_line=False)
        windows = [smoothing._effective_window(call.args[1], len(line.coords))
                   for call in calls.call_args_list]
        self.assertEqual(windows, [4, 3, 2])


if __name__ == '__main__':
    unittest.main()
