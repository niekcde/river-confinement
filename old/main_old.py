#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# """
# Created on Tue Apr 30 10:36:11 2024
# @author: niekcollotdescury
# """

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
from datetime import datetime as dt
from tqdm import tqdm 
import sys

import geopandas as gpd
import pandas as pd
import numpy as np
from multiprocessing import Pool

# %% Import custom modules
# 
os.chdir(directory + 'python/py_code')

from remove_mannual_add import remove_man_add
from reach_def import create_ud_id, create_new_segments
from inflection_points import inflection_points
from width_ratio import width_ratio
from dem import find_dem
from get_orthogonals import get_orthogonals
from cycle_identification import find_connected_side
from support import create_dir,  get_local_utm_projection

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






#%% main
import psutil
import gc
def check_memory():
    process = psutil.Process(os.getpid())
    print(f"Memory usage: {process.memory_info().rss / (1024 * 1024):.2f} MB")


# from memory_profiler import profile 
# @profile
def main(file):

        contName   = file[-29:-27] 
        fileNumber = file[-10:-8]
        print(f'Start: {contName}_{fileNumber}. Time: {dt.now()}')
        code_failure = 0
        
        projection       = 'EPSG:3857'
        projection_coord = 'EPSG:4326'
        
        create_new = False
        calc       = True

        ########################################################################
        # create new reach segments
        ########################################################################
        def create_new_reach_files(directory, c, i):
            vector_file = glob.glob(directory + f"input/SWOT_vector/{contName}*{fileNumber}_*shp")[0]
            node_file   = glob.glob(directory + f'input/SWOT_nodes/{contName}*{fileNumber}_*shp')[0]
        
            df_vector       = gpd.read_file(vector_file)
            if df_vector.crs is None:
                df_vector = df_vector.set_crs(projection_coord)

            df_vector = df_vector.to_crs(projection)
            df_vector_whole = df_vector.copy()

            # add local projection to dataframe
            df_vector_whole = get_local_utm_projection(df_vector_whole)
            
            df_vector_whole = df_vector_whole[df_vector_whole['type'] != 6]      # exclude ghostreaches
            df_vector_whole = df_vector_whole[df_vector_whole.reach_len > (500)] # minimum reachLength


            # finds connected reaches and creates cycles
            cycleFile = directory + f'results/cycles/{c}_{fileNumber}'
            df_cycle, dfConnectionNodes  = find_connected_side(df_vector_whole, projection, cycleFile)
            
            
            df_vector_whole = pd.merge(df_vector_whole, df_cycle, how = 'left', on = 'reach_id')
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
            df_new_segments = create_new_segments(df_upstream, df_upstream, df_node, projection)
            df = df_new_segments.copy()

            # width ratio
            df = width_ratio(df, df_node)

            # Save files
            df.to_file(vector_save_file      
                        , driver='ESRI Shapefile', crs = projection)
            df_node.to_file(node_save_file      
                        , driver='ESRI Shapefile', crs = projection)
            
            return df, df_node
        
        try:
            vector_save_file = directory +f'results/new_segments/vector/{contName}_{fileNumber}_new_segments_vector.shp'
            node_save_file   = directory +f'results/new_segments/node/{contName}_{fileNumber}_new_segments_node.shp' 
            if (len(glob.glob(vector_save_file)) == 1) & (len(glob.glob(node_save_file)) == 1) &\
                (create_new == False):
                try:
                    df       = gpd.read_file(vector_save_file)
                    df_node  = gpd.read_file(node_save_file)

                except:
                    df, df_node = create_new_reach_files(directory, contName, fileNumber)

            else:
                df, df_node = create_new_reach_files(directory, contName, fileNumber)

        except Exception as e:
            print(f"Error in opening file {contName} {fileNumber}: {e}", file=sys.stderr)
 


        ########################################################################
        # Calculate metrics
        ########################################################################
        if (calc == True):
            try:

                #%% set all columns to be filled with values to none
                addColumns = ['infP5_0', 'apex5_0', 'apexP5_0', 'ang5_0', 
                              'outSlope', 'innSlope', 'outPR', 'innPR' , 'CDists', 
                              'bendLines', 'bendWidths', 'bendWidthM']
                df[addColumns] = None

                
                #%% Calc values
                ids = df.reach_id.values
                maxSize = 2000
                int(np.ceil(len(ids) / maxSize))
                idsList = [ids[i:i+maxSize] for i in range(0,len(ids),maxSize)]
                dfTotal = df.copy()

                for idsInd, ids in enumerate(idsList):
                    df = dfTotal[dfTotal.reach_id.isin(ids)]

                    
                    # for id in tqdm(ids):
                    for id in ids:

                        # check_memory()


                        reach = df[df.reach_id == id]
                        ############################
                        # inflection point method
                        ############################
                        if (reach.iloc[0].reach_len < (12*4*30)) | (reach.iloc[0].node_mwm < 30):
                            continue

                        try:
                            sin5_0, infP5_0, apex5_0, apexP5_0, ang5_0 = inflection_points(id, df, df_node,
                                                                                                        projection, True,
                                                                                                        s = 5, degree = 0,
                                                                                                        end_points = True)  

                        except:
                            sin5_0, infP5_0, apex5_0, apexP5_0, ang5_0 = np.nan, np.nan, np.nan, np.nan, np.nan

                        
                            code_failure += 1
                        apexP5_0_wkt = [I.wkt for I in apexP5_0.copy()]
                        infP5_0_wkt  = [I.wkt for I in infP5_0.copy()]
                        df.loc[reach.index, 'sin5_0']     = sin5_0
                        df.at[reach.index[0], 'apex5_0']  = apex5_0
                        df.at[reach.index[0], 'ang5_0']   = ang5_0
                        df.at[reach.index[0], 'apexP5_0'] = apexP5_0_wkt
                        df.at[reach.index[0], 'infP5_0']  = infP5_0_wkt
                        # check_memory()
                        # print()
                        ############################
                        # Confinemet method
                        ############################
                        reach = df[df.reach_id == id]
                        # check_memory()
                        cross_slope_distance = 25
                        inf_settings_input   = [5,0]

                        try:
                            increasedBuffer = 2 # Extra buffer size
                            demWidth = reach.iloc[0].node_mwm * increasedBuffer
                            raster = find_dem(reach, directory, projection,
                                                demWidth*cross_slope_distance, False) 
                            
                            apexP = apexP5_0
                            if isinstance(apexP, list) == True:
                                # check_memory()
                                so,si,po,pi,cd, lineSlope, bl, bw, bwm = get_orthogonals(reach,df_node, raster, 
                                                                    projection, demWidth*increasedBuffer, 
                                                                    cross_slope_distance,inf_settings_input, 
                                                                    directory, True, 400)
                                # check_memory()
                            else:
                                so = si = po = pi = cd = lineSlope = bl = bw = np.nan
                        except:
                            so = si = po = pi = cd= lineSlope = bl = bw = np.nan
                            code_failure += 1

                        
                        
                        # check_memory()
                        df.at[reach.index[0], 'outSlope']  = so
                        df.at[reach.index[0], 'innSlope']  = si
                        df.at[reach.index[0], 'outPR']     = po
                        df.at[reach.index[0], 'innPR']     = pi
                        df.at[reach.index[0], 'CDists']    = cd
                        df.at[reach.index[0], 'bendLines'] = bl
                        df.at[reach.index[0], 'bendWidths']= bw
                        df.at[reach.index[0], 'bendWidthM']= bwm
                        df.loc[reach.index, 'lineSlope']   = lineSlope
                        # check_memory()
                        # gc.collect()
                    df.to_csv(directory + f'temp_files/{contName}_{fileNumber}_{idsInd}.csv'    
                        , index=False)

                tempFiles = glob.glob(directory + f'temp_files/{contName}_{fileNumber}*')
                for ti, f in enumerate(tempFiles):
                    dfS = pd.read_csv(f)
                    if ti == 0:
                        dfTotal = dfS
                    else:
                        dfTotal = pd.concat([dfTotal, dfS], ignore_index=True)

                    os.remove(f) # remove temp File
                
                if cross_slope_distance < 10:
                    dfTotal.to_csv(directory + f'results/all/{contName}_{fileNumber}_0{cross_slope_distance}_{inf_settings_input[0]}_{inf_settings_input[1]}.csv'    
                        , index=False)
                else:
                    dfTotal.to_csv(directory + f'results/all/{contName}_{fileNumber}_{cross_slope_distance}_{inf_settings_input[0]}_{inf_settings_input[1]}.csv'    
                        , index=False)
            


                print(f'Finish: {contName}_{fileNumber} with cross slope {cross_slope_distance}, code failures: {code_failure}. Time: {dt.now()}')
            
            except Exception as e:
                print(f"Error in loop for file {contName} {fileNumber}: {e}", file=sys.stderr)
                
            return df


#%%

#%%

files = glob.glob(directory + 'input/SWOT_vector/*.shp')
files = np.sort(files)

F = glob.glob(directory + 'input/SWOT_vector/*.shp')
for f in range(len(F)):
    D = gpd.read_file(F[f])
    D = D[(D['type'] == 1) & (D.lakeflag == 0)]

    cont = F[f][-29:-27]
    cont_reg = F[f][-10:-8]
    D.loc[:, 'File'] = f'{cont}_{cont_reg}'
    
    if f == 0:
        dfInput = D.copy()
    else:
        dfInput = pd.concat([dfInput, D])

dfInput = dfInput.groupby('File', as_index=False).size()
dfInput = dfInput.sort_values('File')
sizeInput = (-dfInput['size'].values).argsort()
files = files[sizeInput]


# check_memory()
# tm1 = dt.now()
# files = [files[2]]
# for f in files:

#     cont = f[-29:-27] 
#     cont_reg = f[-10:-8]
#     print(cont, cont_reg)
#     df = main(f)
#     print(dt.now() - tm1)
# check_memory()
        



#%%
if __name__ == '__main__':
    with Pool(15) as p:
        p.imap(main, files)
        p.close()
        p.join()



#%%

# try at the top
# try around open
# except at the end

# [2806+892+3000::] + 350


######### final position AS_43