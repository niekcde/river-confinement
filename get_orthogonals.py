

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:36:54 2024

@author: niekcollotdescury
"""
# import packages
import numpy as np
import shapely
import geopandas as gpd
import gc

# import partial packages
from datetime import datetime as dt
from shapely.geometry import LineString
from osgeo import gdal

# import custom modules
from line_functions import create_angled_lines, adjust_confinement_line
from extract_slope_along_raster_line import extract_slope_along_raster_line
from calc_functions import extend_apex
from dem import get_raster_vrt
from support import check_memory


def remove_trailing_missing(arr, missing):
    arr  = np.asarray(arr)  # Ensure it's a NumPy array
    arr[arr == missing] = np.nan
    mask = np.isnan(arr)   # Identify missing values

    if not np.any(mask):   # No missing values
        return arr
    elif np.all(mask): # all missing valyes
        return np.nan
    else: # partial missing values
        last_non_missing_idx = np.where(~mask)[0][-1]  
        if np.any(mask[:last_non_missing_idx+1]): # missing not at the end
            return np.nan
        else: # only missing values at the end
            return arr[:last_non_missing_idx + 1]  

def get_slope_values(slope_line, raster, slope_samples, demFillValeu):
                
    slope_line_angleP, slope_line_angleM = create_angled_lines(slope_line, 50)


    D, P     = extract_slope_along_raster_line(raster, 
                                                    slope_line, 
                                                    samples = slope_samples) 

    D_P, P_P = extract_slope_along_raster_line(raster, 
                                                slope_line_angleP, 
                                                samples = slope_samples) 

    D_M, P_M = extract_slope_along_raster_line(raster, 
                                                slope_line_angleM, 
                                                samples = slope_samples) 
    
    # Mean height profile line to mitigate dem measurement errors
    
    
    elevationValues = []
    for elev in [P, P_P, P_M]:

        elev = remove_trailing_missing(elev, demFillValeu)
        if isinstance(elev, np.ndarray):
            elevationValues.append(elev)

    if len(elevationValues) == 0:
        meanP, P, P_P, P_M,D, D_P, D_M = 99999, 99999, 99999, 99999, 99999, 99999, 99999
    else:
        min_length = min(len(arr) for arr in elevationValues)

        # Trim all arrays to the shortest length
        elevationValues = [arr[:min_length] for arr in elevationValues]
        D               = D[:min_length]
        meanP = np.mean(elevationValues, axis = 0)
    
    return (meanP, P, P_P, P_M, D, D_P, D_M)

def get_orthogonals(line,df:'gpd.GeoDataSeries',
                    apex_points, apexO_points, amplitudes, infLines,
                    bendLines, bendWidths,
                    reachCRS:str, DEMprojection,
                    cross_distance:int,
                    demFillValeu, demVRT,
                    directory:str, slope_samples = 100):
    """
    Function that creates lines at approximately right angles with the river 
    and determines the slope of these lines.\n
    input:
    - line: centerline LineString
    - df: Reach Dataframe (DataSeries)
    - Raster: DEM raster clipped to river segment extent
    - apex_points: list of apex points in local projeciton
    - apexO_points: list of apex origin points on the inflection line
    - inf_points: list of inflection points
    - bendLines: list of bendLines
    - bendWidths: list of bendWidths
    - reachCRS: local utm projection
    - DEMprojection: projection of DEM
    - width_input: reach max width --> check for raster size
    - cross_distance: the factor that half the width is multiplied with to determine the length of 
                      the line created at right angels with the river. 
    - demFillValue: fill value for the selected DEM
    - demVRT: virtual dataset for DEM
    - directory: folder directory
    - slope_samples: Number of samples extracted in slope function\n
    
    return:
    - Elevation profile outerbend
    - Elevation profile innerbend
    - distance outbend
    - distance innerbend 
    - linestring of outerbend
    - linestring of innerbend
    - binary for left right. 1 is left, 0 is right
    - centerlineslope\n

    If confining slope line intersects with the centerline the line is cutoff at the centerline. 
    This can create differences in length between the inner and outer bend
    """
    ##############################
    missingVal = 99999


    ##############################
    # open raster
    ##############################
    increasedBuffer = 1.2 # Extra buffer size
    demWidth        = np.max(bendWidths) * increasedBuffer
    # raster,_ = find_dem_FAB(df, demWidth*cross_distance, dfDemBounds, 
    #                 reachCRS, DEMprojection)
    raster = get_raster_vrt(demVRT, df, demWidth*cross_distance, reachCRS, DEMprojection)
    if isinstance(raster, float):
        return  np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    rb = shapely.box(*raster.rio.bounds()).exterior  

    # set minimum starting length for river segment
    if (line.length > 0) & (isinstance(apex_points, float) == False) &\
        (isinstance(infLines, float) == False): #also checked in the main function

        loopSize = len(apex_points)
        elevOut, elevInn  = np.empty(loopSize, dtype=object), np.empty(loopSize, dtype=object)
        distOut, distInn  = np.empty(loopSize, dtype=object), np.empty(loopSize, dtype=object)
        lineOut, lineInn  = np.empty(loopSize, dtype=object), np.empty(loopSize, dtype=object)
        leftRight         = np.empty(loopSize)
        bendHeight        = np.empty(loopSize)
        for i in range(loopSize): # loop over apex points

            AP        = line.interpolate(line.project(apex_points[i])) # project apex point onto the centerline
            APO       = apexO_points[i] # select apex origin point on inflection line

            # select bendWidth and bendLine
            bendWidth = bendWidths[i]
            bendLine  = bendLines[i]

            _, bh = extract_slope_along_raster_line(raster, bendLine, samples = slope_samples) 
            
            

            # check for each confinement line if it intersects with raster box
            while AP.distance(rb) < (bendWidth*cross_distance): 
                demWidth *= 1.2
                raster.close()
                # raster,_ = find_dem_FAB(df, demWidth*cross_distance, dfDemBounds, 
                #                   reachCRS, DEMprojection)
                raster = get_raster_vrt(demVRT, df, demWidth*cross_distance, reachCRS, DEMprojection)
                rb = shapely.box(*raster.rio.bounds()).exterior


            ##############################
            # Create confinement lines
            ##############################

            D, E = extend_apex(AP, APO, infLines[i], bendWidth, amplitudes[i], cross_distance)
            slope_line_out = LineString([AP, D])
            slope_line_inn = LineString([AP, E])    

            # check if value is left or right side
            offsetLine = line.offset_curve(30)
            if slope_line_out.intersects(offsetLine):
                leftRight[i] = 1
            else:
                leftRight[i] = 0
            

            slope_line_out = adjust_confinement_line(slope_line_out, line)
            slope_line_inn = adjust_confinement_line(slope_line_inn, line)
            
            lineOut[i] = slope_line_out
            lineInn[i] = slope_line_inn

            ##############################
            # get confinement elevation values and distances
            ##############################
            (meanPOut,
                 P_out, P_outP, P_outM, 
                 D_out, D_outP, D_outM) = get_slope_values(slope_line_out, raster, 
                                                           slope_samples, demFillValeu)

            (meanPInn,
                 P_inn, P_innP, P_innM, 
                 D_inn, D_innP, D_innM) = get_slope_values(slope_line_inn, raster, 
                                                           slope_samples, demFillValeu)
            



            ##############################
            # prepare output
            ##############################
            if isinstance(meanPOut, int):
                elevOut[i] = missingVal
            else:
                elevOut[i] = list(meanPOut)
            if isinstance(meanPInn, int):
                elevInn[i] = missingVal
            else:
                elevInn[i] = list(meanPInn)
            
            if isinstance(D_out, int):
                distOut[i] = missingVal
            else:
                distOut[i] = list(D_out)

            if isinstance(D_inn, int):
                distInn[i] = missingVal
            else:
                distInn[i] = list(D_inn)
            
            if np.isnan(np.nanmean(bh)):
                bendHeight[i] = missingVal
            else:
                bendHeight[i] = np.nanmean(bh)
        
        ##############################
        # calculate centerline elevation and slope
        ##############################
        lineDist, lineElev = extract_slope_along_raster_line(raster, 
                                                    line, 
                                                    samples = slope_samples) 

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


    else:
        (elevOut, elevInn, distOut, distInn,
         lineOut, lineInn, leftRight,
         lineSlope, bendHeight) = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan



    raster.close()
    del raster, rb
    gc.collect()
    leftRight = np.array(leftRight)
    return (elevOut, elevInn, distOut, distInn,
            lineOut, lineInn, leftRight,
            lineSlope, bendHeight)

