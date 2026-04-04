# imort packages
import rioxarray
import numpy as np
import xarray
import warnings
import geopandas as gpd
import pandas as pd
import xarray as xr
import os
from pathlib import Path

# import partial packages
from rioxarray.merge import merge_arrays
from shapely import box as Box
from osgeo import gdal

# import tqdm
import gc
#%%


def _collect_fabdem_files(fabdem_dir) -> list[str]:
    fabdem_dir = Path(fabdem_dir).expanduser().resolve()
    dem_files = sorted(
        str(path)
        for path in fabdem_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in {".tif", ".tiff"}
    )

    if len(dem_files) == 0:
        raise FileNotFoundError(f"No FABDEM GeoTIFF files found in {fabdem_dir}")

    return dem_files


def build_fabdem_vrt(fabdem_dir, vrt_path, overwrite: bool = False) -> Path:
    vrt_path = Path(vrt_path).expanduser().resolve()
    vrt_path.parent.mkdir(parents=True, exist_ok=True)

    if vrt_path.exists():
        if overwrite is False:
            return vrt_path
        os.remove(vrt_path)

    dem_files = _collect_fabdem_files(fabdem_dir)
    vrt_ds = gdal.BuildVRT(str(vrt_path), dem_files)
    if vrt_ds is None:
        raise RuntimeError(f"Failed to build FABDEM VRT at {vrt_path}")
    vrt_ds.FlushCache()
    vrt_ds = None
    return vrt_path


def build_fabdem_bounds(fabdem_dir, bounds_path, dem_crs: str = 'EPSG:4326', overwrite: bool = False) -> gpd.GeoDataFrame:
    bounds_path = Path(bounds_path).expanduser().resolve()
    bounds_path.parent.mkdir(parents=True, exist_ok=True)

    if bounds_path.exists() and overwrite is False:
        return gpd.read_file(bounds_path)

    dem_files = _collect_fabdem_files(fabdem_dir)
    ids = []
    geometries = []
    for dem_file in dem_files:
        raster = rioxarray.open_rasterio(dem_file, cache=False)
        try:
            raster_crs = str(raster.rio.crs) if raster.rio.crs is not None else dem_crs
            if raster_crs != dem_crs:
                bounds = raster.rio.transform_bounds(dem_crs)
            else:
                bounds = raster.rio.bounds()
        finally:
            raster.close()

        ids.append(dem_file)
        geometries.append(Box(*bounds))

    df_bounds = gpd.GeoDataFrame({'id': ids, 'geometry': geometries}, crs=dem_crs)
    if bounds_path.exists():
        os.remove(bounds_path)
    df_bounds.to_file(bounds_path, driver='GPKG')
    return df_bounds

def find_dem_bounds_FAB(directory=None, demCRS='EPSG:4326', create_new = False,
                        fabdem_dir=None, dem_boundary_file=None):
    if fabdem_dir is None:
        if directory is None:
            raise ValueError("Either 'directory' or 'fabdem_dir' must be provided.")
        fabdem_dir = Path(directory).expanduser().resolve() / 'input' / 'FAB_dem'

    if dem_boundary_file is None:
        if directory is None:
            raise ValueError("Either 'directory' or 'dem_boundary_file' must be provided.")
        dem_boundary_file = Path(directory).expanduser().resolve() / 'input_created' / 'dem' / 'FAB_dem_bounds.gpkg'

    return build_fabdem_bounds(
        fabdem_dir=fabdem_dir,
        bounds_path=dem_boundary_file,
        dem_crs=demCRS,
        overwrite=create_new,
    )

def find_dem_FAB(rowIn:gpd.GeoDataFrame, bufferSize:int, dfDemBounds,
                 localCRS, demCRS:str='EPSG:4326') -> xarray.core.dataarray.DataArray:
    """
    find matching FAB dem tif file\n
    input:
    - rowIn: input geodataframe row (single row)
    - bufferSize: buffer for centerline
    - directory: main directory
    - demCRS: DEM crs --> default EPSG:4326
    output:
    - raster: rioxarray raster in centerline CRS
    """
    # print('start')
    # print(f'\tFind DEM FAB 1: {check_memory()}')
    row = rowIn.copy(deep = True)
    # print(f'\tFind DEM FAB 2: {check_memory()}')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        row.loc[rowIn.index, 'geometry'] = row.geometry.buffer(bufferSize)
    # print(f'\tFind DEM FAB 3: {check_memory()}')
    
    row    = row.to_crs(demCRS)
    bounds = row.geometry.total_bounds
    # print(bufferSize)
    # print(bounds)
    rowBox = Box(*bounds)
    # print(f'\tFind DEM FAB 4: {check_memory()}')
    
    boolInter = rowBox.intersects(dfDemBounds.geometry)
    dfDemRows = dfDemBounds[boolInter]

    # print(f'\tFind DEM FAB 5 {dfDemRows.shape[0]}: {check_memory()}')
    # print(dfDemRows.shape[0])
    if dfDemRows.shape[0] > 0:
        rasters = []
        for i, demRow in dfDemRows.iterrows():    
            f = demRow['id']
            rasterOpen = rioxarray.open_rasterio(f, cache=False)
            rasters.append(rasterOpen)
        if len(rasters)>1:
            # Merge/Mosaic multiple rasters using merge_arrays method of rioxarray
            res1 = rasters[0].rio.resolution()
            rasterReturn = merge_arrays(dataarrays = rasters,
                                        res = (res1[0], res1[0]), crs=rasters[0].rio.crs)
            del res1
        else:
            rasterReturn = rasters[0]


        rasterReturn = rasterReturn.rio.clip_box(*bounds)
        rasterReturn = rasterReturn.rio.reproject(localCRS)
        
        for r in rasters:
            r.close()
            del r
        del rasters
    else:
        print('ERROR NO MATCHING DEM FOUND')
        rasterReturn = np.nan
    
    
    del row, bounds, boolInter
    # print(f'\tFind DEM FAB 6: {check_memory()}')
    gc.collect()

    
    return rasterReturn, dfDemRows


def get_raster_vrt(vrt, dfRIn,bufferSize, localCRS, demCRS):

    dfR = dfRIn.copy(deep = True)
    # print(f'\tFind DEM FAB 2: {check_memory()}')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dfR.loc[dfRIn.index, 'geometry'] = dfR.geometry.buffer(bufferSize)

    dfR    = dfR.to_crs(demCRS)
    bounds = dfR.geometry.total_bounds

    # Define bounding box (xmin, ymax, xmax, ymin) in the VRT's coordinate system
    bounding_box = (bounds[0], bounds[3], bounds[2], bounds[1])  # bounds in DEM crs

    cropped_ds = gdal.Translate('', vrt, projWin=bounding_box, format='VRT', outputType=gdal.GDT_Float32)
    

    reproj_ds = gdal.Warp('', cropped_ds, dstSRS=localCRS, format='VRT')

    # Read raster data as numpy array
    band = reproj_ds.GetRasterBand(1)
    raster_array = band.ReadAsArray()

    # Get geotransform (needed for coordinate mapping)
    gt = reproj_ds.GetGeoTransform()
    xmin, xres, _, ymax, _, yres = gt
    xmax = xmin + (reproj_ds.RasterXSize * xres)
    ymin = ymax + (reproj_ds.RasterYSize * yres)
    # Create coordinate arrays
    x_coords = np.linspace(xmin, xmax, reproj_ds.RasterXSize)
    y_coords = np.linspace(ymax, ymin, reproj_ds.RasterYSize)  # y is reversed
    # Create xarray dataset
    
    xarr = xr.DataArray(raster_array, coords=[y_coords, x_coords], dims=["y", "x"])
    
    cropped_ds.FlushCache()
    reproj_ds.FlushCache()
    band, reproj_ds, cropped_ds, raster_array = None, None, None, None


    del reproj_ds, raster_array,cropped_ds
    del xmin, xres, ymax, yres
    del xmax, ymin, x_coords, y_coords
    gc.collect()
    return xarr
