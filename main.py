#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# """
# Created on Tue Apr 30 10:36:11 2024
# @author: niekcollotdescury
# """
print('start Code')
directory = '/scratch/6256481/'
# directory = '/Users/niekcollotdescury/Desktop/PhD/RiverWidth/'
# directory = '/Users/6256481/Library/CloudStorage/OneDrive-UniversiteitUtrecht/Desktop/RivWidth/'

#%%
# Supress all python warnings --> is necessary due to many unavoidable shapely warnings
import warnings
warnings.filterwarnings("ignore")



#%% Import packages
print('load Packages')
import glob
import os
from datetime import datetime as dt
from tqdm import tqdm 
import sys
import gc
import re

import geopandas as gpd
import pandas as pd
import numpy as np
from multiprocessing import Pool
from shapely.geometry import LineString
from osgeo import gdal

# %% Import custom modules
# 
os.chdir(directory + 'python/py_code')
from inflection_points import inflection_points_curve
from get_orthogonals import get_orthogonals
from support import create_dir, SWORD_stats, smooth_factor,\
                    file_sorting, node_position, adjust_new_segments,\
                    check_memory
from reach_definition import new_reach_definition
from connect_geometries import merge_centerlines
from smoothing import SG_smoothing

# %%create needed directories

create_dir(directory + 'input_created')
create_dir(directory + 'input_created/dem')
create_dir(directory + 'results')
create_dir(directory + 'results/new_segments')
create_dir(directory + 'results/new_segments/node')
create_dir(directory + 'results/new_segments/vector')
create_dir(directory + 'results/all')
create_dir(directory + 'results/centerline')
create_dir(directory + 'results/cycles')



print('startCode')

os.environ["MALLOC_MMAP_MAX_"] = "40960"

#%% main
def main(multiInput):
        file, confFactor = multiInput

        contName   = file[-29:-27] 
        fileNumber = file[-26:-24]
        print(f'Start: {contName}_{fileNumber}. Time: {dt.now()}')
        code_failure = 0
        
        projection       = 'EPSG:3857'
        projection_coord = 'EPSG:4326'
        DEMprojection    = 'EPSG:4326'
        
        vectorFile = glob.glob(directory + f'results/new_segments/vector/{contName}_{fileNumber}_*.gpkg')[0]
        nodeFile   = glob.glob(directory + f'results/new_segments/node/{contName}_{fileNumber}_*.gpkg')[0]
        df     = gpd.read_file(vectorFile)
        dfNode = gpd.read_file(nodeFile)
        print('read file finished')
        # change rch up/dn from str to list
        df = adjust_new_segments(df) # change string values to list values


        ########################################################################
        # Calculate metrics
        ########################################################################
        minLenFactor = 4*12


        #%% Calc values
        ids = df.loc[df['include_flag'] == '0','combined_reach_id'].unique()
        dfSF = pd.read_csv(directory + 'results/smoothingFactor.csv')
        
        
        vrt_file = directory + "input_created/FAB_dem_vrt.vrt"
        demVRT   = gdal.Open(vrt_file)
        
        # array to string replacement
        str_replace = '[()]|list'

        print(f'Start for loop {len(ids)}')
        
        # for id in tqdm(ids):
        for id in ids:
            calculated = '0'
            try:
                dfReach      = df[df['combined_reach_id'] == id].copy()
                dfReachNodes = dfNode[dfNode['reach_id'].isin(dfReach['reach_id'].values)].copy()

                GroupedCRS = dfReach.groupby('localCRS', as_index = False).size()
                reachCRS   = GroupedCRS.loc[GroupedCRS['size'] == GroupedCRS['size'].max()
                                            ,'localCRS'].iloc[0]

                dfReach      = dfReach.to_crs(reachCRS)
                dfReachNodes = dfReachNodes.to_crs(reachCRS)


                combinedLine, _, _ = merge_centerlines(dfReach, df, reachCRS)

                dfReachNodes       = node_position(combinedLine, dfReachNodes)


                maxWidth = dfReach['combined_reach_max_width'].iloc[0]
                width    = dfReach['combined_reach_width'].iloc[0]
                reachLen = dfReach['combined_reach_len'].iloc[0]
                minLen       = minLenFactor * reachLen

                if np.isnan(maxWidth):
                    maxWidth = dfReachNodes['max_width'].mean()

                ###########################
                # Width Ratio
                ###########################
                widthRatio = dfReachNodes['width'] / dfReachNodes['max_width']
                df.loc[dfReach.index, 'widthRatio'] = np.nanmean(widthRatio)

                ###########################
                # Smoothing
                ###########################
                factorRow       = abs(dfSF['combined_reach_width'] - width).argsort()
                smoothingFactor = dfSF.loc[factorRow, 'smoothFactor'].iloc[0]
                
                smoothing_window = smoothingFactor*int(maxWidth)             # Smoothing window based on mean node max width
                    

                # Smooth the initial line and change dataframe to Series
                combinedLine = SG_smoothing(combinedLine, smoothing_window, maxWidth)
                if len(combinedLine.coords) < 3:
                    combinedLine = combinedLine.segmentize(1) 
                
                df.loc[dfReach.index, 'combined_reach_crs'] = reachCRS           

            except:
                calculated += '1'
                print(id, calculated)

            ############################
            # inflection point method
            ############################

            try:
                (sin,bendSin, infP, infPTotal,
                    apex, apexP, apexPO,
                    ang,
                    bendLines, infLines,
                    bendWidths, bendMaxWidths, bendDistOut, bendLen) = inflection_points_curve(combinedLine, 
                                                dfReach, dfReachNodes) 
                
                # Single values
                bendSin           = np.array2string(bendSin, separator = ', ')
                ang               = np.array2string(ang, separator = ', ')
                bendDistOut       = np.array2string(bendDistOut, separator = ', ')
                bendLen           = np.array2string(bendLen, separator = ', ')
                
                apex_str          = np.array2string(apex, separator = ', ')
                bendWidths_str    = np.array2string(bendWidths, separator = ', ')
                bendMaxWidths_str = np.array2string(bendMaxWidths, separator = ', ')
                
                # geometries
                apexP_wkt         = str([I.wkt for I in apexP.copy()])
                bendLines_wkt     = str([I.wkt for I in bendLines.copy()])
                infLines_wkt      = str([I.wkt for I in infLines.copy()])               
                
            except:
                sin,bendSin, apex, apexP, apexPO, ang = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                infLines, bendLines                   = np.nan, np.nan 
                bendWidths, bendMaxWidths             = np.nan, np.nan 
                bendDistOut, bendLen                  = np.nan, np.nan
                
                apex_str, bendWidths_str, bendMaxWidths_str = np.nan, np.nan, np.nan
                apexP_wkt, infLines_wkt, bendLines_wkt      = np.nan, np.nan, np.nan
                
                code_failure += 1
                calculated   += '2'
            
            # check update method
            df.loc[dfReach.index, 'sin']           = sin                 # Reach Sinuosity 
            df.loc[dfReach.index, 'bendSin']       = bendSin             # Bend Sinuosity
            df.loc[dfReach.index, 'apex']          = apex_str            # Amplitude
            df.loc[dfReach.index, 'apexP']         = apexP_wkt           # Apex Point
            df.loc[dfReach.index, 'ang']           = ang                 # Average bend Curvature
            
            df.loc[dfReach.index, 'infP']          = infLines_wkt        # Inflection Line
            df.loc[dfReach.index, 'bendLines']     = bendLines_wkt       # Bend line
            df.loc[dfReach.index, 'bendWidths']    = bendWidths_str      # Bend Width
            df.loc[dfReach.index, 'bendMaxWidths'] = bendMaxWidths_str   # Bend Max Width
            df.loc[dfReach.index, 'bendDistOut']   = bendDistOut         # Bend Distance to outlet
            df.loc[dfReach.index, 'bendLen']       = bendLen             # Bendline Lenth

            ############################
            # Confinement method
            ############################
            try:
                if isinstance(apexP, np.ndarray) == True:

                    (po,pi,cdo, cdi,
                    lineOut, lineInn, 
                    leftRight, lineSlope, bendH) = get_orthogonals(combinedLine, dfReach,
                                                            apexP, apexPO, apex, 
                                                            infLines, bendLines, 
                                                            bendMaxWidths,
                                                            reachCRS, 'EPSG:4326',
                                                            confFactor, 
                                                            -9999, demVRT,
                                                            directory, 400)

                    po,pi,cdo,cdi    = np.array2string(po, separator=','), np.array2string(pi, separator=','), np.array2string(cdo, separator=','), np.array2string(cdi, separator=',')
                    po, pi = re.sub(str_replace, '', po), re.sub(str_replace, '', pi)
                    cdo, cdi = re.sub(str_replace, '', cdo), re.sub(str_replace, '', cdi)

                    leftRight = np.array2string(leftRight, separator = ', ')
                    bendH     = np.array2string(bendH, separator = ', ')

                else:
                    po = pi = cdo = cdi = lineOut = lineInn = leftRight = lineSlope = bendH = np.nan
                    calculated +='3'

            except:
                po = pi = cdo = cdi = lineOut = lineInn = leftRight = lineSlope = bendH = np.nan
                code_failure += 1
                calculated += '4'
                # print(id, calculated)

            if isinstance(lineOut, np.ndarray):
                lineOut = str([I.wkt if isinstance(I, LineString) else 9999 for I in lineOut])
            if isinstance(lineInn, np.ndarray):
                lineInn = str([I.wkt if isinstance(I, LineString) else 9999 for I in lineInn])

            df.loc[dfReach.index, 'elevOut']    = po
            df.loc[dfReach.index, 'elevInn']    = pi
            df.loc[dfReach.index, 'distOut']    = cdo
            df.loc[dfReach.index, 'distInn']    = cdi
            df.loc[dfReach.index, 'lineOut']    = lineOut
            df.loc[dfReach.index, 'lineInn']    = lineInn
            df.loc[dfReach.index, 'LROrthog']   = leftRight
            df.loc[dfReach.index, 'lineSlope']  = lineSlope
            df.loc[dfReach.index, 'bendHeight'] = bendH
            
            # add value for error or calculation type
            df.loc[dfReach.index, 'calculated']   = calculated

            del dfReach, dfReachNodes
            del apexP,infLines, bendLines


        if confFactor < 10:
            csd = f'0{confFactor}'
            
        else:
            csd = f'{confFactor}'

        fileName = f'{contName}_{fileNumber}_{csd}.csv'    
        df.to_csv(directory +'results/all/' +fileName 
                , index=False)


        print(f'Finish: {contName}_{fileNumber} with cross slope {confFactor}, code failures: {code_failure}. Time: {dt.now()}')
  


#%%
def create_new_reaches_main(c):
    reachFiles = glob.glob(directory + f'input/SWOT_vector/{c}*17*gpkg')
    nodeFiles  = glob.glob(directory + f'input/SWOT_nodes/{c}*17*gpkg')
    # print(reachFiles[0][-25:-23])
    
    print('Create new files')
    for i, fr in enumerate(reachFiles):
        fileName = fr[-25:-23]
        existingFiles = glob.glob(directory + f'results/new_segments/vector/{fileName}*')
        
        for ef in existingFiles:
            print(f'remove: {ef}')
            os.remove(ef)
        print(fileName)
        
        fn = nodeFiles[i]
        
        dfIn     = gpd.read_file(fr)
        dfNodeIn = gpd.read_file(fn)

        new_reach_definition(dfIn,dfNodeIn,4*12, directory, fileName, save = True)

def run_create_new_reaches_main(continents:'list'):
    create_new_stat = True


    print('start multiproces')
    with Pool(len(continents)) as p:
        p.imap(create_new_reaches_main, continents)
        p.close()
        p.join()
    print('Multiproces Done')
    if create_new_stat == True:
        files = np.sort(glob.glob(directory +f'results/new_segments/vector/??_??_*.gpkg'))
        for i, f in enumerate(files):
            D = gpd.read_file(f)
            D['file_cont'] = f[-29:-27]
            D['file_id']   = f[-26:-24]
            if i == 0:
                dfT = D.copy()
            else:
                dfT = pd.concat([dfT, D])
        dfT['file'] = dfT['file_cont'] + '_' + dfT['file_id']
        # dfInc = dfT[dfT['include_flag'] == '0']

        SWORD_stats(dfT, directory)
        smooth_factor(dfT, directory)
        file_sorting(dfT, directory)


# create_new      = True
# if create_new == True:
#     run_create_new_reaches_main(['af', 'oc', 'as', 'sa', 'eu', 'na'])


# # newsegment files
# dfFiles = pd.read_csv(directory + 'results/file_sorting.csv', index_col = 0)
# dfFiles = dfFiles.sort_values('size', ascending = False)
# dfFiles = dfFiles[dfFiles['file'].str.startswith('af')].sort_values('size', ascending = False)

# files = dfFiles['filePath'].values


# files = [files[-1]]
# confFactor = 50
# tm1 = dt.now()
# for f in files:

#     cont = f[-29:-27] 
#     cont_reg = f[-26:-24]
#     print(cont, cont_reg)
#     df = main([f, confFactor])
    # print(dt.now() - tm1)


#%%
continentInput       = sys.argv[1]  # First argument
number_of_processors = int(sys.argv[2])  # Second argument

create_new      = False
if create_new == True:
    run_create_new_reaches_main([continentInput])


dfFiles = pd.read_csv(directory + 'results/file_sorting.csv', index_col = 0)
dfFiles = dfFiles[dfFiles['file'].str.startswith(continentInput)].sort_values('size', ascending = False)
files = dfFiles['filePath'].values
print(files, continentInput, number_of_processors)
print()

confFactor = 50
removeFiles = glob.glob(directory + f'results/all/{continentInput}*_{confFactor}.csv')
for rmf in removeFiles:
    os.remove(rmf)

multiInput = [[f, confFactor] for f in files]
# multiInput = [[files, confFactor]]
# print(multiInput, number_of_processors)
if __name__ == '__main__':
    print('__main__')
    with Pool(number_of_processors) as p:
        p.imap(main, multiInput)
        p.close()
        p.join()


print('code Finished')
