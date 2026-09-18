import math
import unittest

from pipeline.support import _parse_numeric_list_string


class Step2NonfiniteListTests(unittest.TestCase):
    def test_numpy_formatted_nonfinite_values(self):
        values = _parse_numeric_list_string("[1.5, nan, -inf, inf]")
        self.assertEqual(values[0], 1.5)
        self.assertTrue(math.isnan(values[1]))
        self.assertEqual(values[2], -math.inf)
        self.assertEqual(values[3], math.inf)

    def test_geometry_strings_are_unchanged(self):
        self.assertEqual(
            _parse_numeric_list_string("['LINESTRING (0 0, 1 1)']"),
            ['LINESTRING (0 0, 1 1)'],
        )

    def test_other_names_remain_invalid(self):
        with self.assertRaises(ValueError):
            _parse_numeric_list_string("[unexpected]")


if __name__ == "__main__":
    unittest.main()
