import unittest

import geopandas as gpd
import numpy as np
from osgeo import gdal, osr
from shapely.geometry import LineString

from pipeline.dem import get_raster_vrt


class DatelineDemCropTests(unittest.TestCase):
    def test_local_crop_stays_bounded_across_antimeridian(self):
        source = gdal.GetDriverByName("MEM").Create("", 3600, 200, 1, gdal.GDT_Float32)
        source.SetGeoTransform((-180, 0.1, 0, 80, 0, -0.1))
        crs = osr.SpatialReference()
        crs.ImportFromEPSG(4326)
        source.SetProjection(crs.ExportToWkt())
        source.GetRasterBand(1).Fill(42)

        reach = gpd.GeoDataFrame(
            geometry=[LineString([(179.8, 66.8), (-179.8, 66.85)])],
            crs="EPSG:4326",
        ).to_crs("EPSG:32660")
        raster = get_raster_vrt(source, reach, 20_000, "EPSG:32660", "EPSG:4326")

        self.assertLess(raster.size, 100_000)
        self.assertTrue(np.isfinite(raster.values).any())
        self.assertAlmostEqual(float(np.nanmedian(raster.values)), 42)


if __name__ == "__main__":
    unittest.main()
