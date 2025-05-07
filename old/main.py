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
from datetime import datetime


import geopandas as gpd
import pandas as pd
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
from reach_def import create_ud_id, create_new_segments


from inflection_points import inflection_points
from width_ratio import width_ratio

from valley_line import river_dir, inflection_valley
from slope_confinemend import valley_slope_and_confinemend
from remove import remove_dems

from dem import find_dem
from get_orthogonals import get_orthogonals

from cycle_identification import get_cycles
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

#%% regenerate function
import subprocess
import fiona
def regenerate_shx(file_path):
    fixed_file_path = file_path.replace(".shp", "_fixed.shp")
    try:
        subprocess.run(["ogr2ogr", "-skipfailures", "-f", "ESRI Shapefile", fixed_file_path, file_path], check=True)
        return fixed_file_path
    except subprocess.CalledProcessError as e:
        print(f"Failed to regenerate .shx file: {e}")
        return None

def read_shapefile(file_path):
    try:
        gdf = gpd.read_file(file_path)
        return gdf
    except fiona.errors.DriverError as driver_error:
        print(f"DriverError: {driver_error}")
        print("Attempting to regenerate .shx file...")
        fixed_file_path = regenerate_shx(file_path)
        if fixed_file_path:
            try:
                gdf = gpd.read_file(fixed_file_path)
                return gdf
            except Exception as e:
                print(f"Failed to read regenerated shapefile: {e}")
        else:
            print("Could not regenerate the .shx file.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


#%% main
def main(file):

        c = file[-29:-27] 
        i = file[-10:-8]
        code_failure = ''

    # try:
        
        projection    = 'EPSG:3857'
        


        def create_new_reach_files(directory, c, i):
            vector_file = glob.glob(directory + f'input/SWOT_vector/{c}*{i}_*shp')[0]
            node_file   = glob.glob(directory + f'input/SWOT_nodes/{c}*{i}_*shp')[0]
        
            df_vector       = gpd.read_file(vector_file)
            if df_vector.crs is None:
                df_vector = df_vector.set_crs(projection)
            else:
                df_vector = df_vector.to_crs(projection)
            
            df_vector_whole = df_vector.to_crs(projection)

            df_cycle        = get_cycles(df_vector_whole)
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
            df_new_segments = create_new_segments(df_upstream, df_upstream, df_node)

            df = df_new_segments.copy()

            # width ratio
            df = width_ratio(df, df_node)
    
            # Save files
            df.to_file(vector_save_file      
                        , driver='ESRI Shapefile', crs = projection)
            df_node.to_file(node_save_file      
                        , driver='ESRI Shapefile', crs = projection)
            
            return df, df_node
        create_new = False
        calc       = True


        vector_save_file = directory +f'results/new_segments/vector/{c}_{i}_new_segments_vector.shp'
        node_save_file   = directory +f'results/new_segments/node/{c}_{i}_new_segments_node.shp' 
        if (len(glob.glob(vector_save_file)) == 1) & (len(glob.glob(node_save_file)) == 1) &\
            (create_new == False):
            try:
                df       = gpd.read_file(vector_save_file)
                df_node  = gpd.read_file(node_save_file)
            except:
                df, df_node = create_new_reach_files(directory, c, i)
        else:
            df, df_node = create_new_reach_files(directory, c, i)
            
        
        # print('\nCalc new variables: width ratio, inflection sinuosity, regularity, valley sinuosity & valley slope\n')
        if (calc == True):
            #%% set all columns to be filled with values to none
            df['infP5_0'] = None
            df['apex5_0'] = None
            df['apexP5_0'] = None
            df['ang5_0'] = None
            # df['infP10_2'] = None
            # df['infP10_5'] = None
            # df['infP20_2'] = None
            # df['infP20_5'] = None
            
            # df['apex10_2'] = None
            # df['apex10_5'] = None
            # df['apex20_2'] = None
            # df['apex20_5'] = None
            
            # df['apexP10_2'] = None
            # df['apexP10_5'] = None
            # df['apexP20_2'] = None
            # df['apexP20_5'] = None
            
            # df['apexV10_2'] = None
            # df['apexV10_5'] = None
            # df['apexV20_2'] = None
            # df['apexV20_5'] = None
            
            
            df['outSlope'] = None
            df['innSlope'] = None
            df['outPR']    = None
            df['innPR']    = None
            df['CDists']   = None
            
            #%% Calc values
            df = get_local_utm_projection(df)


            ids = df.reach_id.values
            for id in ids:

                reach = df[df.reach_id == id]
            
                if reach.iloc[0].geometry.length < (12*4*30):
                    continue

                try:
                    # inflection_points(ids : 'int', df_vector : 'gpd.GeoPandasDataFrame',
                    #   DFNode : 'gpd.GeopandasDataFrame',
                    #   projection : 'str', plot: 'bool' = True,
                    #   s : 'int' = 5, degree : 'int' = 0, end_points : 'bool' = True):  
                    # return(sin_YE, infP_YE, apex_YE, apexP_YE, wAng_YE)
                    sin5_0,reg5_0, infP5_0, apex5_0, apexP5_0, ang5_0 = inflection_points(id, df, df_node,
                                                                                                  projection, False,
                                                                                                  s = 5, degree = 0,
                                                                                                  end_points = True)  
                    # sin10_2,reg10_2, infP10_2, apex10_2, apexP10_2, apexV10_2 = inflection_points(id, df, df_node,
                    #                                                                               projection, False,
                    #                                                                               s = 5, degree = 0,
                    #                                                                               end_points = True)  
                    # sin10_5,reg10_5, infP10_5, apex10_5, apexP10_5, apexV10_5 = inflection_points(id, df, projection,
                    #                                                                             False, 10, 5, True)
                    # sin20_2,reg20_2, infP20_2, apex20_2, apexP20_2, apexV20_2 = inflection_points(id, df, projection,
                    #                                                                             False, 20, 5, True)
                    # sin20_5 ,reg20_5, infP20_5, apex20_5, apexP20_5, apexV20_5 = inflection_points(id, df, projection,
                    #                                                                             False, 20, 5, True)
                except:
                    sin5_0, reg5_0, infP5_0, apex5_0, apexP5_0, ang5_0 = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                    # sin10_5=reg10_5= infP10_5= apex10_5= apexP10_5= apexV10_5  = np.nan
                    # sin20_2=reg20_2= infP20_2= apex20_2= apexP20_2= apexV20_2  = np.nan
                    # sin20_5=reg20_5= infP20_5= apex20_5= apexP20_5= apexV20_5  = np.nan

                    code_failure += 'Inlfection Point failure'
                
                df.loc[reach.index, 'sin5_0']     = sin5_0
                df.loc[reach.index, 'reg5_0']     = reg5_0
                df.at[reach.index[0], 'infP5_0']  = infP5_0
                df.at[reach.index[0], 'apex5_0']  = apex5_0
                df.at[reach.index[0], 'apexP5_0'] = apexP5_0
                df.at[reach.index[0], 'ang5_0']   = ang5_0

                # df.loc[reach.index, 'sin10_2']    = sin10_2
                # df.loc[reach.index, 'sin10_5']    = sin10_5
                # df.loc[reach.index, 'sinReg10_2'] = reg10_2
                # df.loc[reach.index, 'sinReg10_5'] = reg10_5
                # df.at[reach.index[0], 'infP10_2']   = infP10_2
                # df.at[reach.index[0], 'infP10_5']   = infP10_5
                # df.at[reach.index[0], 'apex10_2']   = apex10_2
                # df.at[reach.index[0], 'apex10_5']   = apex10_5
                # df.at[reach.index[0], 'apexP10_2']  = apexP10_2
                # df.at[reach.index[0], 'apexP10_5']  = apexP10_5
                # df.at[reach.index[0], 'apexV10_2']  = apexV10_2
                # df.at[reach.index[0], 'apexV10_5']  = apexV10_5
                
                # df.loc[reach.index, 'sin20_2']    = sin20_2
                # df.loc[reach.index, 'sin20_5']    = sin20_5
                # df.loc[reach.index, 'sinReg20_2'] = reg20_2
                # df.loc[reach.index, 'sinReg20_5'] = reg20_5
                # df.at[reach.index[0], 'infP20_2']   = infP20_2
                # df.at[reach.index[0], 'infP20_5']   = infP20_5
                # df.at[reach.index[0], 'apex20_2']   = apex20_2
                # df.at[reach.index[0], 'apex20_5']   = apex20_5
                # df.at[reach.index[0], 'apexP20_2']  = apexP20_2
                # df.at[reach.index[0], 'apexP20_5']  = apexP20_5
                # df.at[reach.index[0], 'apexV20_2']  = apexV20_2
                # df.at[reach.index[0], 'apexV20_5']  = apexV20_5

                ########
                # Confinemet Inflection point method
                reach = df[df.reach_id == id]

                cross_slope_distance = 5
                inf_settings_input   = [5,0]
        
                try:
          
                    
                    # CHECKKKKKK
                    raster = find_dem(reach, directory, projection,
                                        reach.iloc[0].node_mwm*cross_slope_distance, False)
                    apexP = apexP5_0
                    # if inf_settings_input == [10,2]:
                    #     apexP = apexP10_2
                    # if inf_settings_input == [10,5]:
                    #     apexP = apexP10_5
                    # if inf_settings_input == [20,2]:
                    #     apexP = apexP20_2
                    # if inf_settings_input == [20,5]:
                    #     apexP = apexP20_5
                    
        
                    if isinstance(apexP, list) == True:
                    
                        so,si,po,pi,cd = get_orthogonals(reach,df_node, raster, 
                                                         projection,reach.iloc[0].max_width, 
                                                         cross_slope_distance,inf_settings_input, 
                                                         directory, False, 100)
                    else:
                        so = si = po = pi = cd = np.nan

                except:
                    so = si = po = pi = cd = np.nan

                    code_failure += '\nConfiment code failure'
                
                df.at[reach.index[0], 'outSlope']  = so
                df.at[reach.index[0], 'innSlope']  = si
                df.at[reach.index[0], 'outPR']     = po
                df.at[reach.index[0], 'innPR']     = pi
                df.at[reach.index[0], 'CDists']    = cd
                
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
        
            # print('\nStart Valley slope calc\n')
            # df['vc'] = df['vc20_2']
            # print(df['vc'])
            # df = valley_slope_and_confinemend(df, df_node, projection, projection_ee, False, directory)
            # remove_dems(df, directory)
            
        # 
            # print('start save files')
            
            # save all geometries to csv 
            if cross_slope_distance < 10:
                df.to_csv(directory + f'results/all/{c}_{i}_{cross_slope_distance}_{inf_settings_input[0]}_{inf_settings_input[1]}.csv'    
                      , index=False)
            else:
                df.to_csv(directory + f'results/all/{c}_{i}_{cross_slope_distance}_{inf_settings_input[0]}_{inf_settings_input[1]}.csv'    
                      , index=False)
            #save file without valley_center
            # df_nv = df.copy()
            # df_nv = df_nv.drop('vc', axis = 1)
            # df_nv.to_file(directory +f'results/centerline/{c}_{i}_centerline.shp'      , driver='ESRI Shapefile', crs = projection)
        

    # except:
    #     time = datetime.now().strftime('%Y-%m-%d %H:%M')
    #     file = open(f'/scratch/6256481/results/all/{c}_{i}_{time}.txt', 'w')
    #     file.write( f"FAILED\n{code_failure}")
    #     file.close()

#%%

#%%
files = glob.glob(directory + 'input/SWOT_vector/*.shp')
files.sort()



files_created = glob.glob(directory + 'results/all/*_10_20_2*')
files_created.sort()

f_input   = [f'{f[-29:-27]}_{f[-10:-8]}' for f in files]
f_created = [f[-17:-12] for f in files_created]


f_not_created = []
for i in range(len(f_input)):

    if f_input[i] not in f_created:
        f_not_created.append(i)

files_not_created = [files[i] for i in f_not_created]
files = files_not_created

# files = [files[f_input.index('sa_62')]]
# files = [files[5]]

# print(len(files))
# for f in files:

#     cont = f[-29:-27] 
#     cont_reg = f[-10:-8]
#     print(cont, cont_reg)
    # main(f)

#%%
# if __name__ == '__main__':
#     with Pool(1) as p:
#         p.map(main, files)
