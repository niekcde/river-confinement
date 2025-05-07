#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# module import
from smoothing import SG_smoothing
from calc_functions import minmax
from line_functions import azimuth_coords, getExtrapoledLine, remove_self_intersect, create_bend_line, get_bend_width
from line_functions import get_points_along_linestring, curvature
from support import change_crs_of_list, get_ids

# package import
import shapely
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



# import of package functions
from shapely.geometry import Point, LineString, MultiPolygon
from scipy.spatial.distance import directed_hausdorff
from datetime import datetime as dt




def inflection_points(ids : 'int', df_vector : 'gpd.GeoPandasDataFrame',
                      DFNode : 'gpd.GeopandasDataFrame',
                      projection : 'str', plot: 'bool' = True,
                      s : 'int' = 5, degree : 'int' = 0, end_points : 'bool' = True):  
    """
    Input:
        ids: row reach ID
        df_vector: vector Dataframe
        DFNode: node dataframe
        projection: Selected projection (standard mercator global)
        plot: boolean value to plot or not, default = True
        s: smoothing window, default = 10
        degree: minimal angle difference for inflection, default = 0
        end_points: include end_points of linestring as inflection points, default = False
    
    return:
        Sinuosity (sin_YE): single value for the sinuosity based on length / inflection length
        Inflection points (infP_YE): inflection points
        Apex values(apex_YE): Apex values per bend
        Apex point per bend (apexP_YE): the location of apex point per bend
        Angle of line (ang_): Averaged difference between points per bend
    """
    #########    #########    #########    #########    #########    #########    #########
    #########    #########    #########    #########    #########    #########    #########
    # start van lijn maakt uit voor de locatie van inflectie punten??
        # Altijd upstream eerste coordinateeeeeeeeee
        # aanpassing moet gemaakt worden in de reach definitie code
    #########    #########    #########    #########    #########    #########    #########
    #########    #########    #########    #########    #########    #########    #########
    dfNode = DFNode.copy()
    
    r           = df_vector[df_vector.reach_id == ids] # select dataframe row
    r           = r.to_crs(r.iloc[0].localCRS)

    vector_orig = r.iloc[0].geometry
    r.loc[:,'geometry'] = remove_self_intersect(r.iloc[0].geometry)

    ###########################
    # Smoothing
    ###########################
    smoothing_window = s*int(r.iloc[0].node_mwm)             # Smoothing window based on mean node max width
    length_ratio     = r.iloc[0].reach_len/r.iloc[0].min_len # Ratio reach length and minimum required length
    
    if length_ratio < 1.0: # If length ratio is below 1 the smoothing window is multiplied by the ratio
        smoothing_window *= length_ratio
        

    # Smooth the initial line and change dataframe to Series
    r.loc[:,'geometry'] = SG_smoothing(r.iloc[0].geometry, smoothing_window, False)
    if len(r.loc[:,'geometry'].iloc[0].coords) < 3:
        r.loc[:,'geometry'] = r.loc[:,'geometry'].iloc[0].segmentize(500) 

    ###########################
    # Select row and id combinations and corresponding Node Rows
    ###########################
    row        = r.iloc[0]
    vector     = row.geometry  # Select geometry
    vec_coords = vector.coords # Get coordinates in the line segment

    rowIDS = get_ids(row)
    if isinstance(rowIDS, list) == False:
        rowIDS = [rowIDS]
    dfNode      = dfNode[dfNode.reach_id.isin(rowIDS)]
    dfNodeL = dfNode.copy()
    dfNodeL = dfNode.to_crs(row.localCRS)

    ########################
    # Calculate azimutal values and cross products for specified points along line
    ########################
    
    # empty lists to be used in loop
    crosses          = []
    fwd              = []
    # loop over all the coordinates minus two. The minus two is to be able to find coordinates 
    # before and after the selected point
    vec_coords = []
    steps = np.linspace(0.005,0.995, int((100/0.5)+1))
    for i in range(len(steps)):
        C = vector.interpolate(steps[i], normalized = True).coords[0]
        vec_coords.append(C)
    
    for c in range(len(vec_coords) -2):

        # Select the target coordinate pair and select point before and after
        target_min  = vec_coords[c]
        target      = vec_coords[c+1]
        target_plus = vec_coords[c+2]
        
        # Calculate direction between 
        fwd_azimuth_minus = azimuth_coords(target_min[0], target_min[1],target[0], target[1])
        fwd_azimuth_plus  = azimuth_coords(target[0], target[1],target_plus[0], target_plus[1])
        fwd.append((fwd_azimuth_minus + fwd_azimuth_plus) / 2)

        # Calculate cross product between points
        v1 = np.array(target) - np.array(target_min)
        v2 = np.array(target_plus) - np.array(target)
        crosses.append((v1[0] * v2[1]) - (v1[1] * v2[0]))
        
    ########################
    # Create dataframe with cross products and coordinate id's
    ########################
    new_df = pd.DataFrame({'id':np.arange(0,len(vec_coords)-2), 'cross':crosses, 'fwd':fwd})

    # get the coordinates for the points 
    new_df[['Xm', 'Ym']] = vec_coords[0:-2]
    new_df[['X', 'Y']]   = vec_coords[1:-1]
    new_df[['Xp', 'Yp']] = vec_coords[2::]
    
    new_df.loc[:,'change']  = np.sign(new_df.cross).diff().ne(0) # check rows where cross product changes sign
    
    new_df.loc[:,'change_pX'] = (new_df.Xm + new_df.X) / 2 # find the point between two points
    new_df.loc[:,'change_pY'] = (new_df.Ym + new_df.Y) / 2 # find the point between two points
    ########################
    # Add first and final coordinate to "changes" data frame for YES dataframe and 
    # remove initial change value for NO dataframe
    ########################
    def adjust_end_points(df):
        df_NO = df.copy()
        df_YE = df.copy()
        
        df_YE.loc[df_YE.index[0], 'change_pX'] = df_YE.loc[df_YE.index[0], 'Xm']
        df_YE.loc[df_YE.index[0], 'change_pY'] = df_YE.loc[df_YE.index[0], 'Ym']

        end_true = df_YE.iloc[-1].change
        if end_true == True:
            df_YE.loc[df_YE.index[-1] + 1, :] = df_YE.loc[df_YE.index[-1], :]
            df_YE.loc[df_YE.index[-1], 'change'] = True
            df_YE.loc[df_YE.index[-1], 'id'] += 1

            df_YE.loc[df_YE.index[-1], 'fwd'] = 1e3
            df_YE.loc[df_YE.index[-1], 'change_pX'] = df_YE.loc[df_YE.index[-1], 'Xp']
            df_YE.loc[df_YE.index[-1], 'change_pY'] = df_YE.loc[df_YE.index[-1], 'Yp']
        else:
            df_YE.loc[df_YE.index[-1], 'change'] = True
            df_YE.loc[df_YE.index[-1], 'change_pX'] = df_YE.loc[df_YE.index[-1], 'Xp']
            df_YE.loc[df_YE.index[-1], 'change_pY'] = df_YE.loc[df_YE.index[-1], 'Yp']
            df_YE.loc[df_YE.index[-1], 'fwd'] = 1e3

    
        df_NO.loc[df_NO.index[0], 'change'] = False
        return df_NO, df_YE
    
    ########################
    # Calculate inflection points from dataframe containing changes in crossproducts
    ########################
    def calc_inflections(df,projection, line,dfNode, width):
        
        new_df = df.copy()
        # Get only values where there is a change in cross-product 
        # --> if end points is true these are changed to true 
        inflections = new_df[new_df.change == True]
        inflections = inflections.copy() # get copy of dataframe to remove chaining errors

        # check the difference in angle with the previous inflection point
        inflections.loc[:,'fwd_change'] = abs(inflections.fwd.diff().values)

        # set the first value to a high value --> no calculation possible
        inflections.loc[inflections.index[0], 'fwd_change'] = 100

        # remove inflection points where the difference in angle with the previous 
        # inflection point is too small
        inflections = inflections[inflections.fwd_change > degree]
        
        gdf_inf     = gpd.GeoDataFrame(
                                        inflections, 
                                        geometry=gpd.points_from_xy(inflections.change_pX, 
                                                                    inflections.change_pY), 
                                        crs=projection)
        
        # check distance between inflection points is X times the river width
        i1   = gdf_inf.index[0]
        inds = [i1]
        for i2 in gdf_inf.index[1::]:
            
            P1 = gdf_inf.loc[i1].geometry          
            P2 = gdf_inf.loc[i2].geometry
            
            infLine   = LineString([P1, P2])
            
            bendLine  = create_bend_line(infLine, line)
            bendWidth, bandWidthMax = get_bend_width(bendLine, dfNode, plot = False)
            bendWidth = bandWidthMax


            # Compute directed Hausdorff distance from bendLine to inflection line
            apexDist = shapely.hausdorff_distance(bendLine, infLine, 0.5)
            
            
            pdist = abs(line.project(P1) - line.project(P2))
            # if inflection point distance more than 4*bendwidth AND apex distance larger than bendwidth add the inflection point
            if (pdist > (bendWidth*4)) & (apexDist > bendWidth):
                inds.append(i2)
                i1 = i2
            
            if (i2 == gdf_inf.index[-1]) & (i1 != i2) & (len(inds) > 1):
                if (pdist < (bendWidth*4)):
                    inds[-1] = i2
                else:
                    inds.append(i2)
            elif (i2 == gdf_inf.index[-1]) & (len(inds) == 1):
                inds.append(i2)

        gdf_inf_filtered = gdf_inf.copy()
        gdf_inf_filtered = gdf_inf_filtered.loc[inds]

        return gdf_inf_filtered    

    ########################
    # Function for sinuosity, ang & Apex
    ########################
    def sin_ang_apex(GDF_INF, vector):
        gdf_inf = GDF_INF.copy()
        if gdf_inf.shape[0] > 1: # only continue if 
        
            ########################
            # Calculate Sinuosity
            ########################
            infLine = LineString(gdf_inf.geometry)
            sinuosity = vector.length / infLine.length


            ########################
            # Angle & Apex 
            ########################
            
            apex            = []  # Apex distance per bend
            apex_points     = []  # apex point per bend
            ang             = [] # List of average angle per bend
            bendLength      = []
            
            for i in range(1, gdf_inf.shape[0]): # Loop over all inflection points

                l = LineString([gdf_inf.iloc[i-1].geometry, gdf_inf.iloc[i].geometry]) 
                
                bendLine = create_bend_line(l, vector)
                ##########
                # Determine apex points
                ##########
                bendLine_coords     = get_points_along_linestring(bendLine, 20)
                infLine_coords      = get_points_along_linestring(infLine, 20)
                apexDist, idxBendLine, idxInfLine = directed_hausdorff(bendLine_coords,
                                                                  infLine_coords)
                apex.append(apexDist)

                apex_points.append(Point(bendLine_coords[idxBendLine]))

                ##########
                # Calculate average angle between points
                ##########

                ang.append(curvature(bendLine))

                bendLength.append(bendLine.length)

            # change all inflection point dataframe in to list
            infP = list(gdf_inf.geometry.values) 
        else:
            # If no inflection points located --> straight line
            sinuosity      = np.nan 
            infP           = np.nan
            apex           = np.nan 
            apex_points    = np.nan 
            ang            = np.nan
            
        return sinuosity, infP, apex, apex_points, ang, bendLength
    
    ########################
    # Run functions
    ########################

    dfNoEnd, dfYeEnd = adjust_end_points(new_df)

    inflectionsYE    = calc_inflections(dfYeEnd, projection,vector,dfNodeL, row.node_mwm)

    sin_YE, infP_YE, apex_YE, apexP_YE, ang_YE, bendL_YE = sin_ang_apex(inflectionsYE, vector)

    
    ########################
    # Determine weightAveraged angle value
    ########################
    wAng_YE = np.average(ang_YE, weights = bendL_YE)
    wAng_YE = ang_YE
    ########################
    # Plot
    ########################
    # if plot == True:
    #     f, ax = plt.subplots(nrows=1, ncols=1, figsize = [10,6], sharey=True)
    #     t1 = f"Length: {int(vector.length)}, Width: {int(row.node_mwm)}, Window: {s}"
    #     t2 = f"Sine: {np.round(sin_YE, 3)}"
    #     t3 = f"Angle: {np.round(np.mean(wAng_YE), 5)}, Apex: {int(np.mean(apex_YE))}"
    #     f.suptitle(f'ID: {ids}\n{t1}\n{t2}, {t3}', fontsize=12)


    #     ax.plot(*vector.xy, color = 'orange', zorder = 0)
    #     ax.plot(*vector_orig.xy, color = 'chocolate', zorder = 0)

    #     def plot_inf_ap(ax, infP, apexP):
    #         if isinstance(infP, list):
    #             inf_line = LineString(infP)
    #             ax.plot(*inf_line.xy, c = 'g', zorder = 10)

    #             for p in infP:
    #                 ax.scatter(*p.xy, c = 'r', zorder= 20, s = 50)

    #             for ap in apexP:
    #                 ax.scatter(*ap.xy, marker = '^', c = 'b', s = 50, zorder = 20)
    #     plot_inf_ap(ax, infP_YE, apexP_YE)

    #     ax.set_aspect('equal', adjustable='box')
    #     ax.axis('off')



    #     plt.tight_layout()
    #     plt.show()

    ########################
    # Change inflection and apex points to global crs
    ########################
    infP_YE  = change_crs_of_list(infP_YE , row.localCRS, projection)
    apexP_YE = change_crs_of_list(apexP_YE, row.localCRS, projection)

    return(sin_YE, infP_YE, apex_YE, apexP_YE, ang_YE)
