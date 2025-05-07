from support import confinement_factor, str_to_list, expand_dataframe, str_to_list_comb
from calc_functions import confinement_values
from confinement_margin import confinement_margin, get_raster_side
from connect_geometries import merge_centerlines

import numpy as np
import os
import pandas as pd
import shapely
import geopandas as gpd
import xarray as xr

from scipy.spatial import KDTree
from tqdm import tqdm 
from osgeo import gdal 
from glob import glob
gdal.UseExceptions()

def ER_slope_margin_values(dfInc,demVRT,directory, cf = [50,10]):
    # singleF = glob(directory + f'results/single_values/*{confSize}.csv')
    # dfE     = open_single_values(singleF)
    # dfInc = dfE[(dfE['include_flag'] == '0') & (dfE['calculated'] == '0000')]
    # dfInc = confinement_factor(dfInc, cf[0], cf[1])
    print(f'run_confinement_values - calc_confinement_values - ER_slope_margin_values: Start')
    
    # Find confinement factor value
    dfCF = pd.read_csv(directory + 'results/confinement_factor.csv')
    # Using KDTree for efficient nearest neighbor search
    tree = KDTree(dfCF[['bendWidths']].values)
    _, idx = tree.query(dfInc[['bendWidths']].values)
    # Getting the closest values
    dfInc['conFactor'] = dfCF.iloc[idx]['conFactor'].values

    
    # change to combined reach id!!!!!!!! so that local crs is correct
    # reachIDS = dfInc['reach_id'].values
    reachIDS = dfInc['combined_reach_id'].astype(int).unique()
    crGeoms = np.empty(len(reachIDS), dtype = object)
    # for i, rid in tqdm(enumerate(reachIDS), total = len(reachIDS)):
    for i, rid in enumerate(reachIDS):

        dfReach     = dfInc[dfInc['combined_reach_id'] == rid].copy()
        # GroupedCRS = dfReach.groupby('localCRS', as_index = False).size()
        # reachCRS   = GroupedCRS.loc[GroupedCRS['size'] == GroupedCRS['size'].max()
        #                                     ,'localCRS'].iloc[0]
        try:
            dfSingle = dfReach.groupby('reach_id').agg({'reach_order':'first', 'geometry':'first', 'dn_connected_reach':'first'})
            dfSingle['geometry'] = shapely.from_wkt(dfSingle['geometry'])
            L, _, _= merge_centerlines(dfSingle, _, _, False)
            crGeoms[i] = L
        except:
            crGeoms[i] = shapely.geometry.LineString()


        for i, r in dfReach.iterrows():
            dfBend = gpd.GeoDataFrame({'geometry':shapely.from_wkt(r['bendLines'])}
                                      , index = [0], crs = r['combined_reach_crs'])
            dfBend = dfBend.to_crs('EPSG:4326')
            bendGeom = dfBend['geometry'].iloc[0]

            (vwo, vwi, cph, cmh, so, si, 
                ERo, ERi) = confinement_values(r['elevOut'], r['elevInn'],
                                                r['distOut'], r['distInn'], 
                                                r['bendWidths'],r['bendMaxWidths'],
                                                r['conFactor'] * 2)
            
            dfInc.loc[i, ['valley_width_out', 'valley_width_inn', 
                    'cp_height', 'cm_height', 
                    'slope_out', 'slope_inn',
                    'ER_out', 'ER_inn', 'bendGeom'
                    ]] = (vwo, vwi, cph, cmh, so, si, ERo, ERi, bendGeom)



    dfInc['ER_left']  = np.where(dfInc['LROrthog'] == 1, dfInc['ER_out'], dfInc['ER_inn'])
    dfInc['ER_right'] = np.where(dfInc['LROrthog'] == 1, dfInc['ER_inn'], dfInc['ER_out'])

    dfInc['slope_left']  = np.where(dfInc['LROrthog'] == 1, dfInc['slope_out'], dfInc['slope_inn'])
    dfInc['slope_right'] = np.where(dfInc['LROrthog'] == 1, dfInc['slope_inn'], dfInc['slope_out'])

    dfInc['ER_max']    = dfInc[['ER_out', 'ER_inn']].max(axis = 1)
    dfInc['ER_min']    = dfInc[['ER_out', 'ER_inn']].max(axis = 1)

    dfInc['slope_max'] = dfInc[['slope_out', 'slope_inn']].max(axis = 1)
    dfInc['slope_min'] = dfInc[['slope_out', 'slope_inn']].max(axis = 1)

    return dfInc, crGeoms

def calc_confinement_values(df,fileName, directory, returnDataframe, open_seperate = False):
    print(f'run_confinement_values - calc_confinement_values: {fileName}')
    vrt_file = directory + "input_created/FAB_dem_vrt.vrt"
    demVRT   = gdal.Open(vrt_file)

    if open_seperate == True:
        listCols = ['distOut', 'distInn', 'elevInn', 'elevOut']
        df = str_to_list(df, listCols, listCols)

    dfE, crGeoms = ER_slope_margin_values(df, demVRT, directory, [50,10])

    ########################################
    # Save reach averaged gpkg
    ########################################
    dfEG = dfE.groupby('combined_reach_id').agg({
    'ER_out':'mean'    ,'ER_inn':'mean',
    'ER_left':'mean'   ,'ER_right':'mean',
    'ER_max':'mean'    ,'ER_min':'mean',
    'slope_out':'mean' ,'slope_inn': 'mean',
    'slope_left':'mean','slope_right': 'mean',
    'slope_max':'mean' ,'slope_min':'mean',
    'cp_height':'mean',
    'catchment_position':'first'}).copy()


    gdfEG = gpd.GeoDataFrame(dfEG, geometry = crGeoms, crs = 'EPSG:4326')
    gdfEG.to_file(directory + f'results/reach_averaged/{fileName}.gpkg', driver = 'GPKG')
    
    ########################################
    # Save reduced bend gpkg
    ########################################
    saveCols = ['reach_id', 'combined_reach_id',
                'apex', 'bendWidths', 'bendMaxWidths',
                'ang', 'sin', 'bendSin','cp_height',
                'ER_out', 'ER_inn','ER_left', 'ER_right',
                'slope_out', 'slope_inn', 'slope_left', 'slope_right',
                'bendDistOut']
    dfE_bend = gpd.GeoDataFrame(dfE[saveCols].copy(), geometry = dfE['bendGeom'], crs= 'EPSG:4326')
    dfE_bend.to_file(directory + f'results/single_values/{fileName}_conf.gpkg', driver = 'GPKG')


    dfE = dfE.drop(['elevOut', 'elevInn', 'distOut', 
                    'distInn', 'geometry', 'bendGeom'], axis = 1)
    dfE['apex'] = dfE['apex'].astype('float64')
    # save nc File
    dfIncX = dfE.to_xarray()
    ncFile = directory + f'results/single_values/{fileName}_conf.nc'
    if os.path.exists(ncFile):
        os.remove(ncFile)
    
    dfIncX.to_netcdf(ncFile)

    if returnDataframe == True:
        return dfE

def concat_nc_conf_files(directory, cross):
    dsList = []
    for c in ['af', 'as', 'sa', 'na', 'oc', 'eu']:
        print(c)
        files = glob(directory + f'results/single_values/{c}*{cross}*conf.nc')
        for f in tqdm(files):
            dsTemp = xr.open_dataset(f)
            dsTemp['file'] = ('index', [c] * dsTemp.sizes['index'])
            dsTemp = dsTemp.drop_vars(['infP', 'bendLines', 'apexP', 'lineInn', 'lineOut'])
            dsList.append(dsTemp)
        
    ds = xr.concat(dsList, dim='index')
    fName = directory + f'results/single_values/global_{cross}_conf.nc'
    if os.path.exists(fName):
        os.remove(fName)
    ds.to_netcdf(fName)

def concat_reachAveraged(directory, cross):
    # for c in ['sa', 'as', 'na', 'eu', 'oc']:
    for c in ['af']:
        print(c)
        files = reachAveraged   = np.sort(glob(directory + f'results/reach_averaged/{c}*{cross}*.gpkg'))
        for i,f in enumerate(files):
            fid = f[-13:-8]
            dfT = gpd.read_file(f)
            dfT['file'] = fid
            if i == 0:
                dfW = dfT
            else:
                dfW = pd.concat([dfW, dfT])
        dfW.to_file(directory + f'results/reach_averaged/{c}_{cross}_reachAveraged_conf.gpkg')


def create_apex_val_dataframe(file, directory = '/scratch/6256481/', returnDataframe = False):

    df = pd.read_csv(file, dtype = {'include_flag': str, 'calculated': str})
    print(df.shape[0])
    dfInc = df[(df['include_flag'] == '0') & (df['calculated'] == '0')]
    dropCols = ['x', 'y', 'swot_obs', 'edit_flag', 'trib_flag', 'wse_var', 'width_var']
    dfInc = dfInc.drop(dropCols, axis = 1)
    print(dfInc.shape[0])

    listStrCols = ['LROrthog', 'bendMaxWidths', 'bendWidths', 'ang', 'bendSin', 'apex', 'bendDistOut', 'bendLen', 'bendHeight'] # new run needed
    # listStrCols = ['LROrthog', 'bendMaxWidths', 'bendWidths', 'ang', 'bendSin', 'apex', 'bendDistOut']
    listGeomCols = ['apexP', 'lineInn', 'lineOut', 'bendLines']
    listNestCols = ['distOut', 'distInn', 'elevInn', 'elevOut']

    listCols = listStrCols + listGeomCols + listNestCols

    dfInc = str_to_list_comb(dfInc, listCols, listNestCols)
    dfE   = expand_dataframe(dfInc.copy())
    
    print('Start saving:', file[-12::])
    dfE.to_csv(directory + f'results/single_values/{file[-12::]}')

    # calc_confinement_values(dfE, file[-12:-4], directory, returnDataframe)