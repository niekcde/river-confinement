import tempfile
import unittest
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.geometry import LineString

from pipeline.sample_step2_profiles import _iter_orthogonal_groups


class OrthogonalBatchReaderTests(unittest.TestCase):
    def test_reaches_split_between_batches_remain_whole_and_ordered(self):
        ids = [1, 1, 2, 2, 2, 3, 3]
        frame = gpd.GeoDataFrame(
            {
                "combined_reach_id": ids,
                "line_type": ["centerline"] * len(ids),
                "geometry": [LineString([(i, 0), (i + 1, 0)]) for i in range(len(ids))],
            },
            crs="EPSG:4326",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "orthogonals.gpkg"
            frame.to_file(path, driver="GPKG")
            actual = list(_iter_orthogonal_groups(path, batch_size=3))

        expected = list(frame.groupby("combined_reach_id", sort=False))
        self.assertEqual([reach_id for reach_id, _ in actual], [reach_id for reach_id, _ in expected])
        for (_, actual_group), (_, expected_group) in zip(actual, expected):
            pd.testing.assert_frame_equal(
                actual_group.reset_index(drop=True), expected_group.reset_index(drop=True)
            )


if __name__ == "__main__":
    unittest.main()
