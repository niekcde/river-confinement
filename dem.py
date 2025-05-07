
# imort packages
import glob
import rioxarray
import numpy as np
import xarray
import warnings
import geopandas as gpd
import pandas as pd
import xarray as xr
import os

# import partial packages
from rioxarray.merge import merge_arrays
from shapely.geometry import Polygon
from shapely import box as Box
from osgeo import gdal

# import tqdm
from datetime import datetime as dt
from support import check_memory
import gc
#%%

def df_add_row(df, folder, name, B):

    newRow = {'Folder':folder,'Name': name, 
              'minx':B[0], 'miny':B[1], 
              'maxx':B[2], 'maxy':B[3]}
    df = pd.concat([df, pd.DataFrame([newRow])], ignore_index=True)
    return df


def find_dem_bounds(directory, projection):
    dem_boundary_file = directory + 'input_created/MERIT_dem_bounds.shp'

    if len(glob.glob(dem_boundary_file)) == 1:
        gdfFindDemBounds = gpd.read_file(dem_boundary_file)
    else:
        dem_folders = glob.glob(directory + 'input/MERIT_dem/*0')
        df = pd.DataFrame(columns = ['Folder','Name', 'minx', 'miny', 'maxx', 'maxy'])

        # for folder in tqdm.tqdm(dem_folders): # loop over all the folders for MERIT
        j = 1
        for i, folder in enumerate(dem_folders):
            
            dem_files = glob.glob(f'{folder}/*.tif')

            for j, f in enumerate(dem_files): # loop over all the tiffs in a folder                

                # open raster and reproject to mercator projection
                raster = rioxarray.open_rasterio(f)    
                if raster.rio.crs != 'EPSG:4326':
                    bounds = raster.transform_bounds("EPSG:4326")
                else:
                    bounds = raster.rio.bounds() 
            
                newRow = {'Folder':folder[-7::],'Name': f[-15:-8], 
                            'minx':bounds[0], 'miny':bounds[1], 
                            'maxx':bounds[2], 'maxy':bounds[3]}
                                
                if (i == 0) & (j == 0):
                    df = pd.DataFrame([newRow])
                else:
                    
                    df = pd.concat([df, pd.DataFrame([newRow])], ignore_index=True)
  
        # create box from the min and max coordinates
        df['geometry'] = df.apply(lambda row: Polygon([
            (row['maxx'], row['miny']),
            (row['maxx'], row['maxy']),
            (row['minx'], row['maxy']),
            (row['minx'], row['miny']),
            (row['maxx'], row['miny'])
            ]), axis=1)
    
        # Step 3: Convert the DataFrame to a GeoDataFrame
        gdf = gpd.GeoDataFrame(df, geometry='geometry')
        
        # Set the Coordinate Reference System (CRS) to WGS84 (EPSG:4326)
        gdf.set_crs(epsg=4326, inplace=True)
        
        gdf.to_file(dem_boundary_file, crs = projection)
    

    return gdfFindDemBounds


def find_dem(row_in:gpd.GeoDataFrame, directory:str, projection:str,buffer_size:float, plot:bool):
    '''Find and clip DEM raster for reach
    input: 
    - input row: input dataframe row
    - directory: main directory for code
    - projection: desired end projection
    - buffer_size: size of the buffer around the input geometry
    - plot: boolean for plotting\n
    
    return:
    - raster
    '''
    row_in = row_in.copy(deep = True)
    row_in.loc[row_in.index, 'geometry'] = row_in.buffer(buffer_size)

    row = row_in.copy()
    row = row_in.to_crs('EPSG:4326')
    row = row.iloc[0]

    # include buffer zone
    df_dem = find_dem_bounds(directory, projection)
    

    intersections = row.geometry.intersects(df_dem.geometry)
    
    df_dem_row = df_dem[intersections]
    if df_dem_row.shape[0] == 1:
        f = directory + f'input/MERIT_dem/dem_tif_{df_dem_row.iloc[0].Folder}/{df_dem_row.iloc[0].Name}_dem.tif'
        raster = rioxarray.open_rasterio(f) 
    elif df_dem_row.shape[0] > 1:
        rasters = []
        for i in range(df_dem_row.shape[0]):    

            f = directory + f'input/MERIT_dem/dem_tif_{df_dem_row.iloc[i].Folder}/{df_dem_row.iloc[i].Name}_dem.tif'
            rasters.append(rioxarray.open_rasterio(f))

        # Merge/Mosaic multiple rasters using merge_arrays method of rioxarray
        res1 = rasters[0].rio.resolution()
        raster = merge_arrays(dataarrays = rasters,
                                  res = (res1[0], res1[0]), crs=rasters[0].rio.crs)

    else:
        print('ERROR NO MATCHING DEM FOUND')
    
    bounds = row.geometry.bounds
    
    raster = raster.rio.clip_box(
            minx=bounds[0],
            miny=bounds[1],
            maxx=bounds[2],
            maxy=bounds[3]
            )

    

    return raster

def find_dem_bounds_FAB(directory, demCRS, create_new = False):
    dem_boundary_file = directory + 'input_created/FAB_dem_bounds.gpkg'

    if (os.path.isfile(dem_boundary_file) & (create_new == False)):
        dfBounds = gpd.read_file(dem_boundary_file)
    else:
        demFiles = glob.glob(directory + 'input/FAB_dem/*.tif')
        name, geom = [], []
        for f in demFiles:
            raster =  rioxarray.open_rasterio(f) 
            if raster.rio.crs != demCRS:
                bounds = raster.transform_bounds(demCRS)
            else:
                bounds = raster.rio.bounds() 
            name.append(f)
            geom.append(Box(*bounds))
        dfBounds = gpd.GeoDataFrame({'id':name, 'geometry':geom}, crs = demCRS)
        dfBounds.to_file(dem_boundary_file, driver = 'GPKG')
    return dfBounds


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
