

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:36:54 2024

@author: niekcollotdescury
"""
import psutil
import os
def check_memory():
    process = psutil.Process(os.getpid())
    print(f"Memory usage: {process.memory_info().rss / (1024 * 1024):.2f} MB")


import numpy as np
import matplotlib.pyplot as plt
import shapely
import geopandas as gpd
import rioxarray
import gc

from datetime import datetime as dt
from shapely.geometry import LineString, Point

# from line_functions import remove_man_add_section, getExtrapoledLine, multiline_check, 


from line_functions import multiline_check, getExtrapoledLine, create_angled_lines, get_bend_width, create_bend_line
from extract_slope_along_raster_line import extract_slope_along_raster_line
from smoothing import SG_smoothing

from support import get_ids
from calc_functions import linear_fit, slope_with_intercept, calculate_point_rightsided_triangle
from dem import find_dem
##############################################################################
# - check if orthog line intersects centerline
##############################################################################
# projection       = 'EPSG:3857'
# projection_coord = 'EPSG:4326'
# import glob
# from support import get_local_utm_projection
# from cycle_identification import find_connected_side
# import pandas as pd
# from remove_mannual_add import  remove_man_add
# from width_ratio import width_ratio
# from reach_def import create_ud_id, create_new_segments
# c = 'oc'
# i = 55
# directory = '/scratch/6256481/'

# create_new = False
# vector_save_file = directory +f'results/new_segments/vector/{c}_{i}_new_segments_vector.shp'
# node_save_file   = directory +f'results/new_segments/node/{c}_{i}_new_segments_node.shp' 
# def create_new_reach_files(directory, c, i):
#             vector_file = glob.glob(directory + f'input/SWOT_vector/{c}*{i}_*shp')[0]
#             node_file   = glob.glob(directory + f'input/SWOT_nodes/{c}*{i}_*shp')[0]
        
#             df_vector       = gpd.read_file(vector_file)
#             if df_vector.crs is None:
#                 df_vector = df_vector.set_crs(projection_coord)

#             df_vector = df_vector.to_crs(projection)
#             df_vector_whole = df_vector.copy()

#             # add local projection to dataframe
#             df_vector_whole = get_local_utm_projection(df_vector_whole)
            
#             df_vector_whole = df_vector_whole[df_vector_whole['type'] != 6]      # exclude ghostreaches
#             df_vector_whole = df_vector_whole[df_vector_whole.reach_len > (500)] # minimum reachLength


#             # finds connected reaches and creates cycles
#             cycleFile = directory + f'results/cycles/{c}_{i}_cycle.shp'
#             df_cycle, dfConnectionNodes  = find_connected_side(df_vector_whole, projection, cycleFile)
            
            
#             df_vector_whole = pd.merge(df_vector_whole, df_cycle, how = 'left', on = 'reach_id')
#             df_vector       = df_vector_whole[(df_vector_whole['type'] == 1) & (df_vector_whole['lakeflag'] == 0)]
            
#             df_node = gpd.read_file(node_file)
#             if df_node.crs is None:
#                 df_node = df_node.set_crs(projection)
#             else:
#                 df_node = df_node.to_crs(projection)
            
            
#             df_node = df_node[df_node.reach_id.isin(df_vector.reach_id)]
            

#             # remove parts of segments with manually added nodes
#             df_vector, df_node_nm = remove_man_add(df_vector,df_node)

#             # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                 # CHANGE TO CORRECT LENGTH CHECK!!!!!!!!!

#             # create new segments from upstream and downstream nodes
#             df_upstream     = create_ud_id(df_vector,df_vector_whole,df_node, df_node_nm)
#             df_new_segments = create_new_segments(df_upstream, df_upstream, df_node, projection)

#             df = df_new_segments.copy()

#             # width ratio
#             df = width_ratio(df, df_node)
    
#             # Save files
#             df.to_file(vector_save_file      
#                         , driver='ESRI Shapefile', crs = projection)
#             df_node.to_file(node_save_file      
#                         , driver='ESRI Shapefile', crs = projection)
            
#             return df, df_node
# # def get_orthogonals(line,width, id,raster, projection, num_segments,cross_distance, plot = False):
# if (len(glob.glob(vector_save_file)) == 1) & (len(glob.glob(node_save_file)) == 1) &\
#     (create_new == False):
#     try:
#         df       = gpd.read_file(vector_save_file)
#         df_node  = gpd.read_file(node_save_file)
#     except:
#         df, df_node = create_new_reach_files(directory, c, i)
# else:
#     df, df_node = create_new_reach_files(directory, c, i)
# Inflection point method





from memory_profiler import profile 
# @profile(precision = 10)
def get_orthogonals(df_row_input:'gpd.GeoDataSeries', nodes_all:gpd.GeoDataFrame,raster, projection:str, width_input:int,
                    cross_distance:int, inf_settings:int, directory:str, plot = False, slope_samples = 100):
    """
    Function that creates lines at approximately right angles with the river 
    and determines the slope of these lines.\n
    input:
    - df_row: Dataframe Row (DataSeries)
    - nodes_all: All nodes in dataframe
    - Raster: DEM raster clipped to river segment extent
    - Projection: the projection used (mercator global)
    - width_input: reach max width --> check for raster size
    - cross_distance: the factor that half the width is multiplied with to determine the length of 
                      the line created at right angels with the river. 
    - inf_settings: the settings used for the inflection point run --> 
                    determines the smoothing applied to the centerline
    - directory: folder directory
    - plot: if True use the build in plot otherwise leave out plot
    - slope_samples: Number of samples extracted in slope function\n
    
    return:
    - Slope values for the outer slope (one value per apex point)
    - Slope values for the inner slope (one value per apex point)
    - Elevation profile outerbend
    - Elevation profile innerbend
    - Distances for inner and outer profile\n

    Values that are not calculated due to a intersection with the centerline are given the value of 9999
    """

    # check_memory()
    df_row = df_row_input.copy()
    df_row = df_row.to_crs(df_row.iloc[0].localCRS)
    row    = df_row.iloc[0]

    if inf_settings == [10,2]:
        apex_points = row.apexP10_2
        inf_points  = row.infP10_2
        smoothing_window = 10
    elif inf_settings == [10,5]:
        apex_points = row.apexP10_5
        inf_points  = row.infP10_5
        smoothing_window = 10
    elif inf_settings == [20,2]:
        apex_points = row.apexP20_2
        inf_points  = row.infP20_2
        smoothing_window = 20
    elif inf_settings == [20,5]:
        apex_points = row.apexP20_5
        inf_points  = row.infP20_5
        smoothing_window = 20
    elif inf_settings == [5,0]:
        apex_points = row.apexP5_0
        inf_points  = row.infP5_0
        smoothing_window = 5
    else:
        print('unknown inflection points settings')
    # check_memory()
    apex_points = [shapely.wkt.loads(A) for A in apex_points]
    inf_points  = [shapely.wkt.loads(A) for A in inf_points]


    smoothing_window = smoothing_window*int(row.node_mwm) # define smoothing window based on mean node max width
    length_ratio     = row.reach_len / row.min_len # Length ratio reach length and minimum required length

    if length_ratio < 1.0: # If length ratio is below 1 the smoothing window is multiplied by the ratio
        smoothing_window *= length_ratio

    
    line = SG_smoothing(row.geometry, smoothing_window, False)
    rid  = row.reach_id
    # check_memory()


    # set minimum starting length for river segment
    if (line.length > 0) & (isinstance(apex_points, float) == False) &\
        (isinstance(inf_points, float) == False): #also checked in the main function


        lineReaches = get_ids(row)
        reachNodes = nodes_all[nodes_all.reach_id.isin(lineReaches)] # all nodes with connected reach_id's
        reachNodes = reachNodes.to_crs(row.localCRS)
        max_dif =[]
        # check_memory()
        slopesOut, slopesInn  = [], []
        elevOut, elevInn      = [], []
        dists                 = []
        bendLines, bendWidths, bendWidthsMax = [], [], []

        dfApex = gpd.GeoDataFrame({'ID':1, 'geometry':apex_points}, crs = projection)
        dfInf  = gpd.GeoDataFrame({'ID':1, 'geometry':inf_points},  crs = projection)
        dfApex = dfApex.to_crs(df_row.crs)
        dfInf  = dfInf.to_crs(df_row.crs )
        apex_points = dfApex.geometry.values
        inf_points  = dfInf.geometry.values
        

        raster = raster.rio.reproject(df_row.crs)
        rb = shapely.box(*raster.rio.bounds()).exterior  

        for i in range(len(apex_points)):
            il = LineString([inf_points[i], inf_points[i + 1]])

            AP        = apex_points[i]
            AP        = line.interpolate(line.project(AP)) 
            AP_origin = il.interpolate(il.project(AP))

            # check_memory()
            bendLine               = create_bend_line(il, line)
            widthWet, widthMax     = get_bend_width(bendLine, reachNodes, plot = False)
            width     = widthMax

            while AP.distance(rb) < (width*cross_distance):

                width_input *= 1.2
                raster.close()

                raster = find_dem(df_row_input, directory, projection, width_input*cross_distance, False)
                raster = raster.rio.reproject(df_row.crs)
                rb = shapely.box(*raster.rio.bounds()).exterior
   
            def create_slope_line(AP1, AP2, width, distance, sign):
                extendDistance = width*distance
                if extendDistance > 30000:
                    extendDistance = 30000
                
                # print(extendDistance, shapely.get_coordinates(AP1)[0], shapely.get_coordinates(AP2)[0])
                extended = getExtrapoledLine(shapely.get_coordinates(AP1)[0], 
                                             shapely.get_coordinates(AP2)[0], 
                                             (extendDistance) + (AP1.distance(AP2)*sign),
                                             1, single = True,
                                             factor_distance = 'Distance')

                # print(list(extended.coords))
                extended = shapely.Point(extended.coords[-1])
                if sign == 1:
                    AP = AP2
                else:

                    AP = AP1

                
                # print(AP, extended)
                slope_line = LineString([AP, extended])


                return slope_line
            
            def inner_outer(slope_line,raster, slope_samples):

                slope_line_angleP, slope_line_angleM = create_angled_lines(slope_line, 200)

            
                S, I, D, P         = extract_slope_along_raster_line(raster, 
                                                                     slope_line, 
                                                                     False, ax = '', 
                                                                     color = 'red', add_reg = False, 
                                                                     samples = slope_samples) 
            
                S_P, I_P, D_P, P_P = extract_slope_along_raster_line(raster, 
                                                                     slope_line_angleP, 
                                                                     False, ax = '', 
                                                                     color = 'red', add_reg = False, 
                                                                     samples = slope_samples) 
            
                S_M, I_M, D_M, P_M = extract_slope_along_raster_line(raster, 
                                                                     slope_line_angleM, 
                                                                     False, ax = '', 
                                                                     color = 'red', add_reg = False, 
                                                                     samples = slope_samples) 
                
            
                meanP          = np.mean([P, P_P, P_M], axis = 0)
                SInterceptMean = slope_with_intercept(D, meanP, meanP[0])
                y_fit          = linear_fit(D, SInterceptMean, meanP[0])

                
                return (SInterceptMean, meanP, y_fit, P, P_P, P_M, D, D_P, D_M,
                       slope_line_angleP, slope_line_angleM)

            def get_slope_values(slope_line, centerline, raster, slope_samples):

                slope_line_inter = LineString(
                        [slope_line.interpolate(0.1,normalized = True),
                         slope_line.interpolate(1  ,normalized = True)])
                # print('0')
                if shapely.intersection_all([slope_line_inter, centerline]).is_empty:
                    (SInterceptMean, meanP, y_fit,
                    P, P_P, P_M, 
                    D, D_P, D_M,
                    slope_line_angleP, slope_line_angleM) = inner_outer(slope_line, raster, slope_samples)
                    
                    # print(1)
                    # print(np.max(meanP))
                    # print(np.min(meanP))

                    if (np.max(meanP) >9e3) | (np.min(meanP) < -100):
                        SInterceptMean = meanP = y_fit        = 99999
                        P = P_P = P_M = D = D_P = D_M         = 99999
                        slope_line_angleP = slope_line_angleM = 99999
                        

                else:
                    # print(2)
                    SInterceptMean = meanP = y_fit        = 99999
                    P = P_P = P_M = D = D_P = D_M         = 99999
                    slope_line_angleP = slope_line_angleM = 99999
                
                return (SInterceptMean, meanP, y_fit, P, P_P, P_M, D, D_P, D_M,
                            slope_line_angleP, slope_line_angleM)
            
            # print(AP_origin, AP,inf_points[i], inf_points[i].distance(AP))
            # check_memory()
            if AP.distance(AP_origin) == 0.0:
                A = shapely.get_coordinates(inf_points[i])[0]
                B = shapely.get_coordinates(AP)[0]
                pR, pL = calculate_point_rightsided_triangle(A[0], A[1], B[0], B[1], AP.distance(inf_points[i]))
                slope_line_out = LineString([AP, Point(pR)])
                slope_line_inn = LineString([AP, Point(pL)])
            else:
                slope_line_out = create_slope_line(AP_origin, AP, width, cross_distance, 1)
                slope_line_inn = create_slope_line(AP, AP_origin, width, cross_distance, 0)

            # check_memory()
            (SInterceptOutMean, meanPOut, y_fit_out,
                 P_out, P_outP, P_outM, 
                 D_out, D_outP, D_outM,
                 slope_line_out_angleP, 
                 slope_line_out_angleM) = get_slope_values(slope_line_out,line, raster, slope_samples)
            # print('INN')
            (SInterceptInnMean, meanPInn, y_fit_inn,
                 P_inn, P_innP, P_innM, 
                 D_inn, D_innP, D_innM,
                 slope_line_inn_angleP, 
                 slope_line_inn_angleM) = get_slope_values(slope_line_inn, line, raster, slope_samples)
            
        
            
            _, _, lineDist, lineElev = extract_slope_along_raster_line(raster, 
                                                        line, 
                                                        False, ax = '', 
                                                        color = 'red', add_reg = False, 
                                                        samples = slope_samples) 
            
            # check_memory()
            # check for missing values in elevation data
            lineDist = np.array(lineDist)
            lineElev[(lineElev < -300 ) | (lineElev == 99999)] = np.nan
            missingElev    = np.isnan(lineElev)
            missingchanges = np.diff(missingElev).sum() # check changes from nan to value

            lineDist = lineDist[missingElev == False]
            lineElev = lineElev[missingElev == False]

            if (missingchanges > 2) | (len(lineElev) == 0): # if there multiple groups of missing data remove slope value
                lineSlope = np.nan
            else:
                lineSlope = np.polyfit(lineDist,lineElev,1)[0]

            # check_memory()
            # if (plot == True) & ( (isinstance(P_inn, float) == False) | (isinstance(P_out, float) == False)) &\
            #    ( (isinstance(meanPOut, int) == False) | (isinstance(meanPInn, int) == False)):
            #     f,(ax1,ax2) = plt.subplots(1,2, figsize=[10,5])
                
    

            #     P = []
            #     titleS2 = ''
            #     titleS3 = ''
            #     if (isinstance(P_out, float) == False) & (isinstance(P_out, int) == False):
            #         # plot 1: elevation profile outerbend
            #         P.extend(P_out)
            #         # Plot height profile vs distance outside bends
            #         ax1.plot(-1*np.array(D_out) , P_out , 'cornflowerblue')
            #         ax1.plot(-1*np.array(D_outP), P_outP, 'cornflowerblue')
            #         ax1.plot(-1*np.array(D_outM), P_outM, 'cornflowerblue')
                
            #         ax1.plot(-1*np.array(D_out), y_fit_out, color = 'g', linestyle = (0, (3, 5, 1, 5)))

            #         S2start = '\nOuterSlope:'
            #         titleS2 = f'{S2start.ljust(20)} {np.round(SInterceptOutMean, 4)}'
                
            #     if (isinstance(P_inn, float) == False) & (isinstance(P_inn, int) == False):
            #         # plot 1: elevation profile innerbend
            #         P.extend(P_inn)

            #         # Plot height profile vs distance inside bends
            #         ax1.plot(D_inn , P_inn , 'navy')
            #         ax1.plot(D_innP, P_innP, 'navy')
            #         ax1.plot(D_innM, P_innM, 'navy')
                
            #         ax1.plot(np.array(D_inn)   , y_fit_inn, color = 'g', linestyle = (0, (3, 5, 1, 5)))

            #         S3start = '\nInner Slope:'
            #         titleS3 = f'{S3start.ljust(20)} {np.round(SInterceptInnMean, 4)}'

            #     # maximal difference in elecation
            #     max_dif.append(np.max(P) -np.min(P))

            #     # create verticale line at centerline position
            #     ax1.vlines(0, np.min(P)-np.min(P)*0.1,np.max(P)+np.max(P)*0.1, color = 'b', linestyle = 'dashed')
                

            #     ax1.set_xlabel('Length (m)')
            #     ax1.set_ylabel('height (m)')
            #     titleS1 = f'Profile with centerline position'
            #     ax1.set_title(f'{titleS1}{titleS2}{titleS3}')


            #     # plot 2: dem overview
            #     ax2.plot(*line.xy, 'lightcoral')
            #     raster[0,:,:].plot.imshow(ax=ax2, cmap = 'RdYlGn_r', 
            #                               vmin=0, 
            #                               vmax=raster.max().values)


            #     if (isinstance(P_out, float) == False) & (isinstance(P_out, int) == False):
            #         ax2.plot(*slope_line_out.xy, 'cornflowerblue')
            #         ax2.plot(*slope_line_out_angleP.xy, 'cornflowerblue',linestyle = '--')
            #         ax2.plot(*slope_line_out_angleM.xy, 'cornflowerblue',linestyle = '--')

             
            #     if (isinstance(P_inn, float) == False) & (isinstance(P_inn, int) == False):
            #         ax2.plot(*slope_line_inn.xy, 'navy')    
            #         ax2.plot(*slope_line_inn_angleP.xy, 'navy',linestyle = '--')
            #         ax2.plot(*slope_line_inn_angleM.xy, 'navy',linestyle = '--')
                
                  
            #     ax2.set_xlabel('x-coordinate (m)')
            #     ax2.set_ylabel('y-coordinate (m)')



            #     T1 = f'{'Elevation raster': ^18}\n{'Width': <14}{int(width)}\n{'Cross:': <12}{int(slope_line_out.length)}\n{'Line Length': <14}{int(line.length)}'
      
            #     ax2.set_title(f'{T1}')
            #     ax2.set_aspect('equal')
                
            #     f.suptitle(rid)
            #     plt.tight_layout()
            #     plt.show()
            
            # check_memory()
            slopesOut.append(SInterceptOutMean)
            slopesInn.append(SInterceptInnMean)

            if isinstance(meanPOut, int):
                elevOut.append(99999)
            else:
                elevOut.append(list(meanPOut))
            if isinstance(meanPInn, int):
                elevInn.append(99999)
            else:
                elevInn.append(list(meanPInn))
            

            if (D_out == 99999):
                if (D_inn != 99999):
                    D = D_inn
                else:
                    D = 99999
            else:
                D = D_out

            dists.append(D) # distance from centerline for elevation measures
            bendLines.append(bendLine.wkt)
            bendWidths.append(widthWet)
            bendWidthsMax.append(widthMax)
            # check_memory()
        # return (slopesOut, slopesInn, elevOut, elevInn, dists, lineSlope, bendLines, bendWidths, bendWidthsMax)
    else:
        slopesOut, slopesInn, elevOut, elevInn, dists, lineSlope, bendLines, bendWidths, bendWidthsMax = \
            np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        
    raster.close()
    del raster
    gc.collect()
    return slopesOut, slopesInn, elevOut, elevInn, dists, lineSlope, bendLines, bendWidths, bendWidthsMax




# from inflection_points import inflection_points
# dfT = df.copy()
# ids = dfT.reach_id.values
# check_memory()
# for id in ids[0:10]:
#     reach = dfT[dfT.reach_id == id].copy()
#     rowIndex = reach.index
#     cross_slope_distance = 60
#     inf_settings_input   = [5,0]
#     check_memory()
#     sin5_0, infP5_0, apex5_0, apexP5_0, ang5_0 = inflection_points(id, dfT, df_node,
#                                                                                         projection, False,
#                                                                                         s = 5, degree = 0,
#                                                                                         end_points = True) 
#     check_memory()
#     demWidth = reach.iloc[0].node_mwm * 2
#     raster = find_dem(reach, directory, projection,
#                         demWidth*cross_slope_distance, False) # *1.2 to create extra buffer
#     check_memory()
#     # raster = raster.rio.reproject(projection)
#     # raster.rio.to_raster(directory + f'results/reach_DEM/{c}_{i}_{id}.tif')
#     # APA = apexP5_0.copy()
#     reach['apexP5_0'] = None
#     reach['infP5_0']  = None

#     reach.at[rowIndex[0], 'apexP5_0'] = [A.wkt for A in apexP5_0]
#     reach.at[rowIndex[0], 'infP5_0']  = [A.wkt for A in infP5_0]

#     # so,si,po,pi,cd, lineSlope, bl, bw, bwm = get_orthogonals(reach,df_node, raster, 
#     #                                     projection, demWidth, 
#     #                                     cross_slope_distance,inf_settings_input, 
#     #                                     directory, False, 400)
#     get_orthogonals(reach,df_node, raster, 
#                                         projection, demWidth, 
#                                         cross_slope_distance,inf_settings_input, 
#                                         directory, False, 400)
#     check_memory()
#     print()

# check_memory()