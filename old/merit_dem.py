import glob

import geopandas as gpd
import pandas as pd
import rioxarray
from rioxarray.merge import merge_arrays
from shapely.geometry import Polygon


def df_add_row(df, folder, name, B):
    newRow = {'Folder': folder, 'Name': name,
              'minx': B[0], 'miny': B[1],
              'maxx': B[2], 'maxy': B[3]}
    df = pd.concat([df, pd.DataFrame([newRow])], ignore_index=True)
    return df


def find_dem_bounds(directory, projection):
    dem_boundary_file = directory + 'input_created/MERIT_dem_bounds.shp'

    if len(glob.glob(dem_boundary_file)) == 1:
        gdfFindDemBounds = gpd.read_file(dem_boundary_file)
    else:
        dem_folders = glob.glob(directory + 'input/MERIT_dem/*0')
        df = pd.DataFrame(columns=['Folder', 'Name', 'minx', 'miny', 'maxx', 'maxy'])

        j = 1
        for i, folder in enumerate(dem_folders):
            dem_files = glob.glob(f'{folder}/*.tif')

            for j, f in enumerate(dem_files):
                raster = rioxarray.open_rasterio(f)
                if raster.rio.crs != 'EPSG:4326':
                    bounds = raster.transform_bounds("EPSG:4326")
                else:
                    bounds = raster.rio.bounds()

                newRow = {'Folder': folder[-7::], 'Name': f[-15:-8],
                          'minx': bounds[0], 'miny': bounds[1],
                          'maxx': bounds[2], 'maxy': bounds[3]}

                if (i == 0) & (j == 0):
                    df = pd.DataFrame([newRow])
                else:
                    df = pd.concat([df, pd.DataFrame([newRow])], ignore_index=True)

        df['geometry'] = df.apply(lambda row: Polygon([
            (row['maxx'], row['miny']),
            (row['maxx'], row['maxy']),
            (row['minx'], row['maxy']),
            (row['minx'], row['miny']),
            (row['maxx'], row['miny'])
        ]), axis=1)

        gdf = gpd.GeoDataFrame(df, geometry='geometry')
        gdf.set_crs(epsg=4326, inplace=True)
        gdf.to_file(dem_boundary_file, crs=projection)

    return gdfFindDemBounds


def find_dem(row_in, directory, projection, buffer_size, plot):
    row_in = row_in.copy(deep=True)
    row_in.loc[row_in.index, 'geometry'] = row_in.buffer(buffer_size)

    row = row_in.copy()
    row = row_in.to_crs('EPSG:4326')
    row = row.iloc[0]

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

        res1 = rasters[0].rio.resolution()
        raster = merge_arrays(
            dataarrays=rasters,
            res=(res1[0], res1[0]),
            crs=rasters[0].rio.crs,
        )
    else:
        print('ERROR NO MATCHING DEM FOUND')

    bounds = row.geometry.bounds
    raster = raster.rio.clip_box(
        minx=bounds[0],
        miny=bounds[1],
        maxx=bounds[2],
        maxy=bounds[3],
    )
    return raster
