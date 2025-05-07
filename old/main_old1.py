#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:36:11 2024

@author: niekcollotdescury
"""


directory = '/scratch/6256481/'
# directory = '/Users/niekcollotdescury/Desktop/PhD/RiverWidth/'
# directory = '/Users/6256481/Library/CloudStorage/OneDrive-UniversiteitUtrecht/Desktop/RivWidth/'

#%%
# Supress all python warnings --> is necessary due to many unavoidable shapely warnings
import warnings
warnings.filterwarnings("ignore")



#%% Import packages
import glob
import os

# import pandas as pd
import geopandas as gpd
import numpy as np
from multiprocessing import Pool
# import scipy as sc

import logging

# import shapely
from shapely.geometry import LineString 


# os.chdir(directory + 'py_code')
os.chdir(directory + 'python/py_code')
# print(os.getcwd())

# from tqdm import tqdm




# %% Import functions
from remove_mannual_add import remove_manual_add, remove_man_add
from smoothing import SG_smoothing
from reach_def import create_ud_id
from reach_def import create_new_segments

from inflection_points import inflection_points
from width_ratio import width_ratio

from valley_line import river_dir, inflection_valley
from slope_confinemend import valley_slope_and_confinemend
from remove import remove_dems



#%%
# create needed directories
def create_dir(path):
    check = os.path.exists(path)
    if check == False:
        os.makedirs(path)
create_dir(directory + 'input_created')
create_dir(directory + 'input_created/dem')
create_dir(directory + 'results')
create_dir(directory + 'results/new_segments')
create_dir(directory + 'results/new_segments/node')
create_dir(directory + 'results/new_segments/vector')
create_dir(directory + 'results/all')
create_dir(directory + 'results/centerline')

#%%

def main(file, DF, DFN):
    
    if (isinstance(DF, gpd.GeoDataFrame)) & isinstance(DFN, gpd.GeoDataFrame):
        all_loaded = True
    else: 
        all_loaded = False
    
    projection    = 'EPSG:3857'
    projection_ee = 'EPSG:4326'
    
    
    openFiles  = True
    create_new = True
    calc = True
    

    c, i = file
    if openFiles == True:

        vector_save_file = directory +f'results/new_segments/vector/{c}_{i}_new_segments_vector.shp'
        node_save_file   = directory +f'results/new_segments/node/{c}_{i}_new_segments_node.shp' 

        if (len(glob.glob(vector_save_file)) == 0) | (len(glob.glob(node_save_file)) == 0) | (create_new == True):

            vector_file = glob.glob(directory + f'input/SWOT_vector/{c}*{i}_*shp')[0]
            node_file   = glob.glob(directory + f'input/SWOT_nodes/{c}*{i}_*shp')[0]
          
            df_vector       = gpd.read_file(vector_file)
            if df_vector.crs is None:
                df_vector = df_vector.set_crs(projection)
            else:
                df_vector = df_vector.to_crs(projection)
            
            df_vector_whole = df_vector.to_crs(projection)
            df_vector       = df_vector_whole[(df_vector_whole['type'] == 1) & (df_vector_whole['lakeflag'] == 0)]
            
            df_node = gpd.read_file(node_file)
            if df_node.crs is None:
                df_node = df_node.set_crs(projection)
            else:
                df_node = df_node.to_crs(projection)
            
            
            df_node = df_node[df_node.reach_id.isin(df_vector.reach_id)]
            

            # remove parts of segments with manually added nodes
            df_vector, df_node_nm = remove_man_add(df_vector,df_node)

            
            # create new segments from upstream and downstream nodes
            df_upstream     = create_ud_id(df_vector,df_vector_whole,df_node, df_node_nm)
            df_new_segments = create_new_segments(df_upstream, df_upstream, df_node)

            df = df_new_segments.copy()

            # width ratio
            df = width_ratio(df, df_node)
    
            # Save files
            df.to_file(vector_save_file      
                        , driver='ESRI Shapefile', crs = projection)
            df_node.to_file(node_save_file      
                        , driver='ESRI Shapefile', crs = projection)

        else:
            if all_loaded == False:
                df       = gpd.read_file(vector_save_file)
                df_node  = gpd.read_file(node_save_file)
            else:
                df = DF.copy()
                df_node = DFN.copy()
    # print('\nCalc new variables: width ratio, inflection sinuosity, regularity, valley sinuosity & valley slope\n')
    if calc == True:
        
    
        # for id in tqdm(df.reach_id.values,
        #           ncols = 50):
        df['infP10_2'] = None
        df['infP10_5'] = None
        df['infP20_2'] = None
        df['infP20_5'] = None
        
        df['apex10_2'] = None
        df['apex10_5'] = None
        df['apex20_2'] = None
        df['apex20_5'] = None
        
        df['apexP10_2'] = None
        df['apexP10_5'] = None
        df['apexP20_2'] = None
        df['apexP20_5'] = None
        
        df['apexV10_2'] = None
        df['apexV10_5'] = None
        df['apexV20_2'] = None
        df['apexV20_5'] = None
        
        ids = df.reach_id.values
        
        for id in ids:
            print(id)
            #try:

            reach = df[df.reach_id == id]
            
            if reach.iloc[0].geometry.length < (12*4*30):
                continue
            
            
            sin10_2,reg10_2, infP10_2, apex10_2, apexP10_2, apexV10_2 = inflection_points(id, df, projection, 
                                                                                          False, 10, 2, True)
            sin10_5,reg10_5, infP10_5, apex10_5, apexP10_5, apexV10_5 = inflection_points(id, df, projection,
                                                                                          False, 10, 4, True)
            sin20_2,reg20_2, infP20_2, apex20_2, apexP20_2, apexV20_2 = inflection_points(id, df, projection,
                                                                                          False, 20, 2, True)
            sin20_5 ,reg20_5, infP20_5, apex20_5, apexP20_5, apexV20_5 = inflection_points(id, df, projection,
                                                                                          False, 20, 4, True)

            df.loc[reach.index, 'sin10_2']    = sin10_2
            df.loc[reach.index, 'sin10_5']    = sin10_5
            df.loc[reach.index, 'sinReg10_2'] = reg10_2
            df.loc[reach.index, 'sinReg10_5'] = reg10_5
            df.at[reach.index[0], 'infP10_2']   = list(infP10_2)
            df.at[reach.index[0], 'infP10_5']   = list(infP10_5)
            df.at[reach.index[0], 'apex10_2']   = apex10_2
            df.at[reach.index[0], 'apex10_5']   = apex10_5
            df.at[reach.index[0], 'apexP10_2']  = apexP10_2
            df.at[reach.index[0], 'apexP10_5']  = apexP10_5
            df.at[reach.index[0], 'apexV10_2']  = apexV10_2
            df.at[reach.index[0], 'apexV10_5']  = apexV10_5
            
            df.loc[reach.index, 'sin20_2']    = sin20_2
            df.loc[reach.index, 'sin20_5']    = sin20_5
            df.loc[reach.index, 'sinReg20_2'] = reg20_2
            df.loc[reach.index, 'sinReg20_5'] = reg20_5
            df.at[reach.index[0], 'infP20_2']   = list(infP20_2)
            df.at[reach.index[0], 'infP20_5']   = list(infP20_5)
            df.at[reach.index[0], 'apex20_2']   = apex20_2
            df.at[reach.index[0], 'apex20_5']   = apex20_5
            df.at[reach.index[0], 'apexP20_2']  = apexP20_2
            df.at[reach.index[0], 'apexP20_5']  = apexP20_5
            df.at[reach.index[0], 'apexV20_2']  = apexV20_2
            df.at[reach.index[0], 'apexV20_5']  = apexV20_5


            # row = reach.iloc[0]
            # vc_20_2 is the correct version
            # vc10_2 = inflection_valley(row, apexP10_2, 50, id,projection, projection_ee, directory, plot = False)
            # vc10_5 = inflection_valley(row, apexP10_5, 50, id,projection, projection_ee, directory, plot = False)
            # vc20_2 = inflection_valley(row, apexP20_2, 50, id,projection, projection_ee, directory, plot = False)
            # vc20_5 = inflection_valley(row, apexP20_5, 50, id,projection, projection_ee, directory, plot = False)
                
            # df.loc[reach.index, 'vc10_2']           = vc10_2
            # df.loc[reach.index, 'vc10_5']           = vc10_5
            # df.loc[reach.index, 'vc20_2']           = vc20_2
            # df.loc[reach.index, 'vc20_5']           = vc20_5
            # print(df.loc[reach.index, 'vc20_2'] , vc20_2)
            # print(type(df))
            
            # if vc10_2.is_empty:
            #     df.loc[reach.index, 'vc_sin10_2'] = np.nan
            # else:
            #     df.loc[reach.index, 'vc_sin10_2'] = row.geometry.length / vc10_2.length
                
            # if vc10_5.is_empty:
            #     df.loc[reach.index, 'vc_sin10_5'] = np.nan
            # else:
            #     df.loc[reach.index, 'vc_sin10_5'] = row.geometry.length / vc10_5.length
            # if vc20_2.is_empty:
            #     df.loc[reach.index, 'vc_sin20_2'] = np.nan
            # else:
            #     df.loc[reach.index, 'vc_sin20_2'] = row.geometry.length / vc20_2.length
            # print(df.loc[reach.index, 'vc_sin20_2'], row.geometry.length / vc20_2.length) 
            # if vc20_5.is_empty:
            #     df.loc[reach.index, 'vc_sin20_5'] = np.nan
            # else:
            #     df.loc[reach.index, 'vc_sin20_5'] = row.geometry.length / vc20_5.length

            # Smooth line before calculating valley line --> 
            # shapely.geometry.LineString.parallel_offset gives multiple linesegments otherwise
            # print(id,2)
            # s = 300
            # reach.loc[:,'geometry'] = SG_smoothing(reach.iloc[0].geometry, s, True)
            # row = reach.iloc[0]
            # node_row = df_node[df_node.reach_id == row.reach_id]
            
            
            # river_direction, off_L, off_R, RLs, buffer,L_length, R_length = river_dir(row,node_row, 50)   

            # init_buffer.append(buffer)
            # init_width.append( node_row.max_width.mean())

            
            # df.loc[reach.index, 'vc']           = river_direction
            # df.loc[reach.index, 'Ll_check']     = L_length
            # df.loc[reach.index, 'Rl_check']     = R_length
            
            # if river_direction.is_empty:
            #     df.loc[reach.index, 'vc_sinus'] = np.nan
            # else:
            #     df.loc[reach.index, 'vc_sinus'] = row.geometry.length / river_direction.length
            #except:
            #    df.loc[reach.index, 'vc20_2']           = LineString()
            #    failed_id.append(id)

            
    #         # print(id,5)
        print('\nStart Valley slope calc\n')
        # df['vc'] = df['vc20_2']
        # print(df['vc'])
        # df = valley_slope_and_confinemend(df, df_node, projection, projection_ee, False, directory)
        # remove_dems(df, directory)
        
    # 
        # print('start save files')
        
        # save all geometries to csv 
        df.to_csv(directory +f'results/all/{c}_{i}.csv'      , index=False)
        
        #save file without valley_center
        # df_nv = df.copy()
        # df_nv = df_nv.drop('vc', axis = 1)
        # df_nv.to_file(directory +f'results/centerline/{c}_{i}_centerline.shp'      , driver='ESRI Shapefile', crs = projection)
    
    if calc == True:
        return df, df_node
    else:
        return df, df_node
#%%

#%%
files = [['na', 71], ['sa', 61], ['as', 45], ['as', 44], ['eu', 21], ['af', 18], ['eu', 23]]
files = glob.glob(directory + 'input/SWOT_vector/*.shp')
# files = [files[0]]

i = 0
import time
for f in files:
    tm = time.time()
    cont = f[-29:-27] 
    cont_reg = f[-10:-8]
    # cont = f[0]
    # cont_reg = f[1]
    
    # print(cont, cont_reg)
    dfV, dfN = main([cont, int(cont_reg)], '', '')
    
    # print(f'total seconds: {time.time() - tm}, {i}')

#%%
# main([cont, int(cont_reg)], dfV, dfN)



# import time
# start_time = time.time()

# print(start_time - time.time())



##### add remove self intersect from reach def and save self intersecting reaches


