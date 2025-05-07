# import packages
import geopandas as gpd
import pandas as pd
import numpy as np
import os
import gc

#import partial packages
from glob import glob
from multiprocessing import Pool

#import custom Modules
from support import check_memory, adjust_new_segments
from dem import find_dem_FAB,find_dem_bounds_FAB

directory = '/scratch/6256481/'

def divide_dataframe_in_equal_parts(df, max_size = 1000):
    new_rows = []
    
    for _, row in df.iterrows():
        value = row['size']
        parts = int(np.ceil(value / max_size))  # Round up division
        step = value // parts  # Divide value equally
        remainder = value % parts  # Handle remainder

        start = 0
        for i in range(parts):
            end = start + step +(remainder if i == (parts-1) else 0)  # Distribute remainder evenly

            new_rows.append({'file':row['file'], 'values': value, 'start': start, 'end': end, 'order':i})
            start = end 

    new_df = pd.DataFrame(new_rows)   
    df = df.merge(new_df, on = 'file', how = 'left')
    return df

def save_raster_reach(dfReach, reachID,reachCRS, cont, bufferSize, dfDemBounds):


    raster, dfDemRows = find_dem_FAB(dfReach, bufferSize, dfDemBounds,reachCRS, 'EPSG:4326')
    if dfDemRows.shape[0] == 0:
        print(dfReach.reach_id.values)
        rasterPath   = ''
        rasterType   = 'None'
        rasterBounds = ''
        rasterCRS    = ''

    elif dfDemRows.shape[0] == 1:
        rasterPath   = dfDemRows.iloc[0]['id']   
        rasterType   = 'single' 
        rasterBounds = str(raster.rio.bounds())
        rasterCRS    = str(raster.rio.crs)
    else:
        rasterPath = directory + f'input_created/FAB_dem_reach/{cont}_{int(reachID)}_DEM.tif'
        raster.rio.to_raster(rasterPath)
        rasterType   = 'multi'
        rasterBounds = str(raster.rio.bounds())
        rasterCRS    = str(raster.rio.crs)

    dfDEM = pd.DataFrame({'combined_reach_id': reachID,'Continent':cont, 
                              'bounds':rasterBounds, 'bounds_crs':rasterCRS,
                              'rasterPath': rasterPath, 'rasterType': rasterType}, index = [0])
    raster.close()
    del raster
    gc.collect()
    return dfDEM

def run_save_raster_reach(df, cont,contOrder, bufferSize, dfDemBounds):
    dfDEM = 0
    ids = df.loc[df['include_flag'] == '0', 'combined_reach_id'].unique()

    for i, rid in enumerate(ids):

        dfReach      = df[df['combined_reach_id'] == rid].copy()
        GroupedCRS = dfReach.groupby('localCRS', as_index = False).size()
        reachCRS   = GroupedCRS.loc[GroupedCRS['size'] == GroupedCRS['size'].max()
                                    ,'localCRS'].iloc[0]
        dfReach      = dfReach.to_crs(reachCRS)


        dfDEMTemp = save_raster_reach(dfReach, rid,reachCRS, cont, 
                                     dfReach['combined_reach_max_width'].iloc[0] * bufferSize,
                                     dfDemBounds)
        if isinstance(dfDEM, int):
            dfDEM = dfDEMTemp
        else:
            dfDEM = pd.concat([dfDEM, dfDEMTemp])


    dfDEM.to_csv(directory + f'input_created/FAB_dem_reach_dataframe/{cont}_{contOrder}_{bufferSize}.csv')

    return dfDEM

def multi_save_raster_reach(multiInput):
    dfDemBounds = find_dem_bounds_FAB(directory, 'EPSG:4326')
    
    filePath, start, end, order = multiInput
    df = gpd.read_file(filePath)
    df = adjust_new_segments(df)
    size = df.loc[df['include_flag'] == '0'].shape[0]
    if end == size:
        dfFilter = df.iloc[start::]
    else:
        dfFilter = df.iloc[start:end]
    print(filePath[-29:-27], order, check_memory())
    run_save_raster_reach(dfFilter, filePath[-29:-27], order, 50, dfDemBounds)
    print(filePath[-29:-27], order, check_memory())
    print()
    del df, dfFilter


if __name__ == '__main__':
    files = np.sort(glob(directory + 'results/new_segments/vector/*'))
    sortFiles = pd.read_csv(directory + 'results/file_sorting.csv')
    sortFiles = sortFiles.sort_values('size')


    newSortFiles = divide_dataframe_in_equal_parts(sortFiles)
    multiInput = newSortFiles[['filePath', 'start', 'end', 'order']].values[0:2]

    files = glob(directory + 'input_created/FAB_dem_reach_dataframe/*')
    [os.remove(f) for f in files]
    print('start Multiprocess', len(multiInput))
    print(multiInput)
    # multi_save_raster_reach(multiInput[0])
    with Pool(1) as p:
        p.imap(multi_save_raster_reach, multiInput)
        p.close()
        # p.terminate()  # Forcefully terminate workers
        p.join()
    
    print('End Multiprocess')
    files = glob(directory + 'input_created/FAB_dem_reach_dataframe/*')
    for i, f in enumerate(files):
        dfTemp = pd.read_csv(f)
        if i == 0:
            dfR = dfTemp
        else:
            dfR = pd.concat([dfR, dfTemp])
    dfR.to_csv(directory + 'input_created/FAB_dem_reach_dataframe/FAB_dem_reach.csv')
    print('code Finished')