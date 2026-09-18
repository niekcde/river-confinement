"""Boundary/category and identity checks for the source/sink analysis."""
import tempfile
import unittest
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import LineString

from analysis.source_sink import join_assignments, overlay, load_vectors, prepare_vectors, CachedRasterSampler, split_dateline


class SourceSinkTests(unittest.TestCase):
    def test_dateline_crossing_is_short_and_ordinary_line_unchanged(self):
        line = LineString([(179.99, 65), (-179.99, 65)])
        fixed = split_dateline(line)
        self.assertEqual(fixed.geom_type, 'MultiLineString')
        self.assertEqual(len(fixed.geoms), 2)
        projected = gpd.GeoSeries([fixed], crs=4326).to_crs(8857)
        self.assertLess(projected.length.iloc[0], 3000)
        normal = LineString([(2, 45), (2.01, 45)])
        self.assertEqual(split_dateline(normal).wkb, normal.wkb)

    def test_cached_sampler_matches_rasterio_at_edges_and_nodata(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'mask_strat_0.tif'
            with rasterio.open(path, 'w', driver='GTiff', width=3, height=1,
                               count=1, dtype='int16', crs='EPSG:8857', nodata=-1,
                               transform=from_origin(0, 100, 100, 100)) as dst:
                dst.write(np.array([[3, -1, 2]], dtype=np.int16), 1)
            coords = np.array([[-1, 50], [0, 100], [99.999, 50], [100, 50],
                               [200, 50], [300, 50], [200, 0], [250, 100.01]])
            with rasterio.open(path) as src:
                np.testing.assert_array_equal(CachedRasterSampler(src).sample(coords),
                                              np.array(list(src.sample(coords))))

    def test_new_fit_mapping_is_explicit(self):
        vectors = pd.DataFrame({'vector_id': [0, 1, 2, 3], 'source_file': ['af_00'] * 4,
                                'bendID': ['0_1', '0_2', '0_3', '0_4'], 'GMM_7_hard': [1, 3, 0, 6]})
        metrics = pd.DataFrame({'vector_id': [0, 1, 2, 3], 'dominant_cat': [3]*4, 'transitions': [0]*4})
        mapping = {1: 'Unconfined', 3: 'Sym partial', 0: 'Asym partial', 6: 'Confined'}
        joined = join_assignments(vectors, metrics, class_map=mapping)
        self.assertEqual(joined.confinement_class.tolist(), list(mapping.values()))

    def test_prepare_and_load_geoparquet_preserves_bend_keys(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            source = gpd.GeoDataFrame({'bendID': ['0_2', '0_1'], 'GMM_7_hard': [0, 2]},
                                      geometry=[LineString([(0, 0), (.01, 0)]),
                                                LineString([(0, 1), (.01, 1)])], crs='EPSG:4326')
            source.to_file(folder / 'eu_00_bend_cluster.gpkg', layer='bends', driver='GPKG')
            output = folder / 'vectors.parquet'
            self.assertEqual(prepare_vectors(folder, '*.gpkg', output), 2)
            loaded = load_vectors(output, geometry=True)
            self.assertEqual(loaded.vector_id.tolist(), [0, 1])
            self.assertEqual(loaded.bendID.tolist(), ['0_2', '0_1'])
            self.assertEqual(loaded.crs.to_epsg(), 8857)
            # Equal Earth is not equidistant: .01 degree is about 958 m here.
            self.assertTrue(loaded.geometry.length.between(950, 970).all())

    def test_three_zones_and_unmatched_bend_keep_original_ids(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'mask_strat_0.tif'
            with rasterio.open(path, 'w', driver='GTiff', width=3, height=1,
                               count=1, dtype='int16', crs='EPSG:8857', nodata=-1,
                               transform=from_origin(0, 100, 100, 100)) as dst:
                dst.write(np.array([[1, 2, 3]], dtype=np.int16), 1)
            vectors = gpd.GeoDataFrame({'vector_id': [17, 99]}, geometry=[
                LineString([(25, 50), (275, 50)]),
                LineString([(400, 50), (450, 50)])], crs='EPSG:8857')
            result, mapping = overlay(vectors, [str(path)])
            self.assertEqual(result.vector_id.tolist(), [17, 99])
            self.assertEqual(mapping.vector_id.tolist(), [17])
            row = result.iloc[0]
            self.assertEqual(row.transitions, 2)
            # Six endpoint samples at 50 m weight: manuscript length is 300, not 250.
            self.assertEqual(row.length_m, 300)
            self.assertEqual(row.dominant_cat, 1)  # first encountered in a three-way tie
            np.testing.assert_allclose(row[['cat_1_pct', 'cat_2_pct', 'cat_3_pct']].astype(float), 1/3)
            self.assertTrue(pd.isna(result.iloc[1].dominant_cat))
            self.assertEqual(result.iloc[1].length_m, 0)

    def test_nodata_gap_legacy_transition_and_coverage(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'mask_strat_0.tif'
            with rasterio.open(path, 'w', driver='GTiff', width=3, height=1,
                               count=1, dtype='int16', crs='EPSG:8857', nodata=-1,
                               transform=from_origin(0, 100, 100, 100)) as dst:
                dst.write(np.array([[3, -1, 2]], dtype=np.int16), 1)
            v = gpd.GeoDataFrame({'vector_id': [5]}, geometry=[LineString([(25, 50), (275, 50)])], crs='EPSG:8857')
            result, _ = overlay(v, [str(path)])
            self.assertEqual(result.iloc[0].length_m, 200)
            self.assertEqual(result.iloc[0].transitions, 1)
            self.assertEqual(result.iloc[0].coverage_pct, 1)  # intentionally historical

    def test_join_reorders_safely_and_rejects_missing_or_duplicate_ids(self):
        vectors = pd.DataFrame({'vector_id': [9, 2], 'source_file': ['as_01', 'eu_00'],
                                'bendID': ['0_1', '0_1'], 'GMM_7_hard': [0, 2]})
        metrics = pd.DataFrame({'vector_id': [2, 9], 'dominant_cat': [1, 3], 'transitions': [0, 0]})
        joined = join_assignments(vectors, metrics).set_index('vector_id')
        self.assertEqual(joined.loc[9, 'zone'], 'Source')
        self.assertEqual(joined.loc[9, 'confinement_class'], 'Unconfined')
        self.assertTrue(joined.figure_eligible.all())
        with self.assertRaises(ValueError):
            join_assignments(vectors, metrics.iloc[:1])
        with self.assertRaises(ValueError):
            join_assignments(vectors, pd.concat([metrics, metrics]))


if __name__ == '__main__':
    unittest.main()
