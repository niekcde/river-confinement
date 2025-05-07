#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# module import
from line_functions import azimuth_coords, create_bend_line, get_bend_width
from line_functions import get_points_along_linestring, curvature, get_bend_dist_out

# package import
import shapely
import geopandas as gpd
import pandas as pd
import numpy as np

# import of package functions
from shapely.geometry import Point, LineString
from scipy.spatial.distance import directed_hausdorff
from datetime import datetime as dt

import matplotlib.pyplot as plt

def inflection_points(vector : 'shapely.geometry.LineString',
                      df : 'gpd.GeoPandasDataFrame',DFNode : 'gpd.GeopandasDataFrame',
                      projection : 'str',
                      degree : 'int' = 0, end_points : 'bool' = True):  
    """
    Input:
        vector: combined reach id

        df: vector Dataframe
        DFNode: node dataframe
        projection: reach projection 
        degree: minimal angle difference for inflection, default = 0
        end_points: include end_points of linestring as inflection points, default = False
    
    return:
        Sinuosity (sin_YE): single value for the sinuosity based on length / inflection length
        Inflection points (infP_YE): inflection points (saved in local CRS)
        Apex values(apex_YE): Apex values per bend
        Apex point per bend (apexP_YE): the location of apex point per bend (saved in local CRS)
        Angle of line (ang_): Averaged difference between points per bend
    """
    dfNode = DFNode.copy()
    
    ###########################
    # select line coordinates
    ###########################
    vec_coords = vector.coords


    ########################
    # Calculate azimutal values and cross products for specified points along line
    ########################
    # empty lists to be used in loop
    crosses          = []
    fwd              = []


    # loop over all the coordinates minus two. The minus two is to be able to find coordinates 
    # before and after the selected point
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
    new_df[['X', 'Y']]      = vec_coords[1:-1]
    new_df.loc[:,'change']  = np.sign(new_df['cross']).diff().ne(0) # check rows where cross product changes sign
    
    ########################
    # Calculate inflection points from dataframe containing changes in crossproducts
    ########################
    def calc_inflections(df,projection, line,dfNode, dfR):
        new_df = df.copy()

        # inflections are points where cross product changes sign
        inflections = new_df[new_df['change'] == True]
        inflections = inflections.copy() # get copy of dataframe to remove chaining errors

        # check the difference in angle with the previous inflection point
        inflections.loc[:,'fwd_change'] = abs(inflections['fwd'].diff().values)

        # set the first value to a high value --> no calculation possible
        inflections.loc[inflections.index[0], 'fwd_change'] = 100
        inflections = inflections[inflections['fwd_change'] > degree]
        # check distance between inflection points is X times the river width 
        gdf_inf     = gpd.GeoDataFrame(inflections['id'], 
                                       geometry=gpd.points_from_xy(inflections['X'],inflections['Y']), 
                                       crs=projection)
        
        infStart = gpd.GeoDataFrame({'id': 0, 'geometry': Point(line.coords[0])}, index = [0])
        infEnd   = gpd.GeoDataFrame({'id': gdf_inf['id'].max() +1,
                                     'geometry': Point(line.coords[-1])}, index = [0])
        gdf_inf = pd.concat([infStart, gdf_inf, infEnd], ignore_index=True)


        i1   = gdf_inf.index[0]
        inds = [i1]
        bendLines, bendWidths, bendMaxWidths = [], [], []
        for i, i2 in enumerate(gdf_inf.index[1::]):
            
            P1 = gdf_inf.loc[i1].geometry          
            P2 = gdf_inf.loc[i2].geometry
            infLine   = LineString([P1, P2])

            bendLine  = create_bend_line(infLine, line)
            bendWidth, bendMaxWidth = get_bend_width(line, bendLine, dfNode, dfR)
            bendWidthCalc = bendMaxWidth 

            # Compute directed Hausdorff distance from bendLine to inflection line
            apexDist = shapely.hausdorff_distance(bendLine, infLine, 0.5)
        

            pdist = abs(line.project(P1) - line.project(P2))
            # if inflection point distance more than 4*bendwidth AND apex distance larger than bendwidth add the inflection point
            if (pdist > (bendWidthCalc*4)) & (apexDist > bendWidthCalc):
                inds.append(i2)
                i1 = i2
                
                bendLines.append(bendLine)
                bendWidths.append(bendWidth)
                bendMaxWidths.append(bendMaxWidth)
                
            if (i2 == gdf_inf.index[-1]) & (i1 != i2) & (len(inds) > 1):
                if (pdist > (bendWidthCalc*4)) | (apexDist > (bendWidthCalc)):
                    inds.append(i2)
                    bendLines.append(bendLine)
                    bendWidths.append(bendWidth)
                    bendMaxWidths.append(bendMaxWidth)
                    
                else:
                    inds[-1] = i2
                    
                    P1 = gdf_inf.loc[inds[-2]].geometry
                    P2 = gdf_inf.loc[i2].geometry

                    infLine   = LineString([P1, P2])
                    bendLine  = create_bend_line(infLine, line)
                    bendWidth, bendMaxWidth = get_bend_width(line, bendLine, dfNode, dfR)

                    bendLines[-1] = bendLine
                    bendWidths[-1] = bendWidth
                    bendMaxWidths[-1] = bendMaxWidth
            
            elif (i2 == gdf_inf.index[-1]) & (len(inds) == 1):
                inds.append(i2)
                bendLines.append(bendLine)
                bendWidths.append(bendWidth)
                bendMaxWidths.append(bendMaxWidth)
            
        gdf_inf_filtered = gdf_inf.copy()
        gdf_inf_filtered = gdf_inf_filtered.loc[inds]
        return gdf_inf_filtered,gdf_inf, bendLines, bendWidths, bendMaxWidths

    ########################
    # Function for sinuosity, ang & Apex
    ########################
    def sin_ang_apex(GDF_INF, vector, bendLines):
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
            apex_points,apexOrigin_points = [],[]  # apex point per bend
            ang             = [] # List of average angle per bend
            bendLength, bendDistOut     = [], []
            
            for i in range(1, gdf_inf.shape[0]): # Loop over all inflection points
                
                inf_section = LineString([gdf_inf.iloc[i-1].geometry, gdf_inf.iloc[i].geometry]) 
                
                bendLine = bendLines[i-1]
                ##########
                # Determine apex points
                ##########
                bendLine_coords     = get_points_along_linestring(bendLine, 20)
                infLine_coords      = get_points_along_linestring(inf_section, 20)
                apexDist, idxBendLine, idxInfLine = directed_hausdorff(bendLine_coords,
                                                                       infLine_coords)
                

                apexPoint = Point(bendLine_coords[idxBendLine])
                apexOriginPoint = Point(infLine_coords[idxInfLine])
                if apexOriginPoint.distance(bendLine) < 10:
                    apexDist = 0

                apex.append(apexDist)
                apex_points.append(apexPoint)
                apexOrigin_points.append(apexOriginPoint)

                ##########
                # Calculate average angle between points
                ##########
                ang.append(curvature(bendLine))
                bendLength.append(bendLine.length)
                bendDistOut.append(get_bend_dist_out(vector, bendLine, dfNode))
                
            # change all inflection point dataframe in to list
            infP = list(gdf_inf.geometry.values) 
        else:
            # If no inflection points located --> straight line
            sinuosity         = np.nan 
            infP              = np.nan
            apex              = np.nan 
            apex_points       = np.nan 
            ang               = np.nan
            apexOrigin_points = np.nan

        return sinuosity, infP, apex, apex_points, apexOrigin_points, ang, bendLength, bendDistOut
    
    ########################
    # Run functions
    ########################
    dfYeEnd = new_df

    (inflectionsYE, dfInfTotal, bendLines, 
     bendWidths, bendMaxWidths)    = calc_inflections(dfYeEnd, projection,vector,dfNode, 
                                        df)

    (sin_YE, 
     infP_YE, apex_YE, apexP_YE, apexPO_YE,
     ang_YE, bendL_YE, bendDO) = sin_ang_apex(inflectionsYE, vector, bendLines)

    ########################
    # Determine weightAveraged angle value
    ########################
    # wAng_YE = np.average(ang_YE, weights = bendL_YE)
    # wAng_YE = ang_YE


    ########################
    # Change inflection and apex points to global crs
    ########################
    # infP_YE  = change_crs_of_list(infP_YE , row.localCRS, projection)
    # apexP_YE = change_crs_of_list(apexP_YE, row.localCRS, projection)
    sin           = sin_YE
    infP          = infP_YE
    infPTotal     = dfInfTotal['geometry'].values
    apex          = np.array(apex_YE)
    apexP         = apexP_YE
    apexPO        = apexPO_YE
    ang           = np.array(ang_YE)
    bendWidths    = np.array(bendWidths)
    bendMaxWidths = np.array(bendMaxWidths)
    bendDO        = np.array(bendDO)
    
    return(sin, infP, infPTotal, apex, apexP, apexPO, ang, bendLines, bendWidths, bendMaxWidths, bendDO)


def get_point_perpindicular_from_point_off_line(line, point_off_line):
    """Calculate Point on straight line that is perpindicular to input point a distance off input line."""
    # Given values
    x1, y1 = line.coords[-1]  # Endpoint A of the line
    x2, y2 = line.coords[0]   # Endpoint B of the line
    x3, y3 = point_off_line.x, point_off_line.y     # Given point P

    # Compute parameter t
    t = ((x3 - x1) * (x2 - x1) + (y3 - y1) * (y2 - y1)) / ((x2 - x1) ** 2 + (y2 - y1) ** 2)

    # Compute the perpendicular point Q(x_q, y_q)
    x_q = x1 + t * (x2 - x1)
    y_q = y1 + t * (y2 - y1)

    return Point([x_q, y_q])

def get_apex_curve(bendLine, infLine, arcS):
    """Get apex point based on a curvature determined value. 
    If line section is classified as straight get apex point as midway point."""
    if arcS == 0:
        apexP     = bendLine.interpolate(0.5, normalized = True)
        apexPO    = Point(apexP)
        amplitude = 0 # pas aan
    else:
        coordArray = np.array(bendLine.coords)
        diffCoords = np.diff(coordArray, axis = 0)

        az = np.rad2deg(np.arctan2(diffCoords[:,1], diffCoords[:,0]))
        az[az<0] = az[az<0]+ 360
        if arcS < 0:
            deltaAz = abs((az[-1]-az[1]) %360)
            adjust = 1
        else:
            deltaAz = (az[0]-az[-1]) % 360
            adjust = -1
        
        if deltaAz <= 180:
            ap = az[0] + ((deltaAz/2)*adjust)
        else:
            ap = az[0] + (90*adjust)

        ap = ap if ap>0 else abs(ap - 360)
        apInd = np.argmin(abs(az - ap))

        apexP  = Point(coordArray[apInd])
        apexPO = shapely.ops.nearest_points(Point(apexP), infLine)[1]
        amplitude = apexP.distance(apexPO)
    return apexP, apexPO, amplitude

def get_apex_distance(bendLine, infLine, arcS):
    """Get apex point based on the maximum minimal distance between two lines (Hausdorf distance). 
    If line section is classified as straight get apex point as midway point.\n
    input:
    - bendLine: LineString of bend
    - infLine: Linestring connection inflection points
    - arcS: Curvature sign of the bend, 0 indicates straight\n
    output:
    - shapely inflection point value of apex
    - shapely inflection point value of apex origin on the inflection line
    - amplitude value"""


    if arcS == 0:
        apexP  = bendLine.interpolate(0.5, normalized = True)
        apexPO = Point(apexP)
        if apexP.distance(infLine) < 1:
            apexPO    = apexP
            amplitude = 0
        else:
            apexPO    = get_point_perpindicular_from_point_off_line(infLine, apexP)
            amplitude = apexP.distance(apexPO)
    else:
        coordArray = np.array(bendLine.coords)

        infLine_coords      = get_points_along_linestring(infLine, 10)
        amplitude, idxBendLine, idxInfLine = directed_hausdorff(coordArray,
                                                                infLine_coords)
        apexP =  Point(coordArray[idxBendLine])
        apexPO = Point(infLine_coords[idxInfLine])
    return apexP, apexPO, amplitude

def find_sign_changes(arr):
    sign_changes = np.where(np.sign(arr[:-1]) != np.sign(arr[1:]))[0]
    return sign_changes

def inflection_points_curve(line:"shapely.LineString", 
                            dfR:"gpd.GeoDataFrame", 
                            dfNodeR:"gpd.GeoDataFrame"):
    """Calculate inflection points based on sign of curvature of a bend. Minimal amplitude requirend of 0.5 bendwidth.\n
    input:
    - line: reach centerline
    - dfR: reach dataframe
    - dfNodeR: reach node dataframe\n
    output:
    - Reach sinuosity
    - Combined Inflection points
    - Initially calculated inlfection points
    - Amplitudes
    - Apex points (list of shapely points)
    - Apex Origin points (list of shapely points)
    - Curvature values per bend
    - Bendlines (list of shapely LineStrings)
    - bendWidths
    - bend Max widths
    - bend distance to outlet
    All output geometries are in local reach crs
    """
    # dt1 = dt.now()
    ##############################
    # Get coordinates and curvature values
    ##############################
    coords = np.array(line.coords)
    curve  = curvature(line, False, False)
    # get list of indices for curvature changes
    curveChanges = find_sign_changes(curve)

    curveChangesVal = curve[curveChanges]
    zeroVals = np.where(curveChangesVal == 0.0)[0]
    curveChanges = np.delete(curveChanges, zeroVals)

    if curveChanges[0] != 0:
        curveChanges = np.insert(curveChanges, 0, 0)
    if curveChanges[-1] != len(curve)-2:
        curveChanges = np.append(curveChanges, len(curve)-1)
    else:
        curveChanges[-1] = len(curve)-1



    ##############################
    # Loop over the initial curvature changes (possible inflection points)
    ##############################
    segments, tp, arcSign, straight = np.empty(len(curveChanges)-1), [],[], False
    p1     = Point(coords[0])
    p1List = [p1]
    curveVal =0
    for i in range(len(curveChanges)-1):
        # dt3 = dt.now()
        # Positions of possible inflection points. Select two points to create a linesegment
        cc1, cc2 = curveChanges[i], curveChanges[i+1] 
        pp1 , p2  = Point(coords[cc1]), Point(coords[cc2])
        infLine = LineString([p1,p2])
        arcLine = create_bend_line(infLine, line)

        arcAmp                  = shapely.hausdorff_distance(infLine, arcLine) # Calculate Arc Amplitde
        bendWidth, bendMaxWidth = get_bend_width(line, arcLine, dfNodeR, dfR) # get arc bendWidth and max bendWidth

        C = curve[cc1]
        C = np.mean(curve[cc1:cc2+1])

        # # get indices for curvature range
        # curveIndex1 = tp[i] if i == 0 else tp[i] + 1
        # curveIndex2 = tp[i+1] +1 if (i +2) < len(tp)  else tp[i+1]
        # print()
        

        # Divide Arc segments in Straight (small amlitude), Positive and Negative
        arcAmpDim = arcAmp / bendMaxWidth
        
        if arcAmpDim <= 0.5:
            segmentSign = 0
        elif C < 0:
            segmentSign = -1
        elif C > 0:
            segmentSign = 1
        else:
            if i == 0:
                segmentSign = np.mean(curvature(arcLine, False, False))
            else:
                segmentSign = segments[i-1]
        # print()
        # print(curve[cc1], C, np.mean(curve[cc1:cc2]))
        segments[i] = segmentSign # add values to list


        # plt.plot(*infLine.xy)
        # plt.plot(*arcLine.xy)
        # plt.axis('equal')
        # plt.show()
        

        # Combine Inflection points. Consequtive inflection points with same sign are combined and 
        # inflection points with same sign with straight section in between are also combined.
        if (i == 0): # always have first coordinate as inflection point
            curveVal = segmentSign
            tp.append(cc1), arcSign.append(segmentSign)
            segments[i] = 0
            # print('\tAdd')
        elif i > 0:
            targetPoint = cc1
            if segmentSign != segments[i-1]:


                p1 = Point(coords[cc2]) # update combined point after sign change in inflection points
                p1List.append(p1)
                
                if segmentSign == 0:
                    straight = True
                    tp.append(targetPoint), arcSign.append(segmentSign)
                    # print('\tstraight add')
                elif (straight == True) & (segmentSign == curveVal):
                    # tp[-1] = targetPoint
                    tp.pop(), arcSign.pop(), p1List.pop()
                    straight = False
                    p1 = p1List[-1]
                    print('\tstraight remove')
                else:
                    tp.append(targetPoint), arcSign.append(segmentSign)
                    # print('\tCurve add')
                    curveVal = segmentSign
                    straight = False
        print(f'{np.round(arcAmpDim, 6)}'.ljust(10), f'{np.round(arcAmp, 5)}'.ljust(9), f'{np.round(C, 5)}'.ljust(9), f'{len(tp)}'.ljust(3), f'{segmentSign}'.ljust(3), curveVal, p1, Point(coords[cc1]),p2)
    if np.where(curveChanges == tp[-1])[0] == (len(curveChanges) -2):
        if (arcAmp / bendMaxWidth) < 0.5:
            tp[-1] = curveChanges[-1]
        else:
            tp.append(curveChanges[-1]) # always add final point as inflection point
    else:
        tp.append(curveChanges[-1]) # always add final point as inflection point
    
    if tp[1] == curveChanges[1]:
        p1, p2 = Point(coords[tp[0]]), Point(coords[tp[1]])

        infLine = LineString([p1,p2])
        arcLine = create_bend_line(infLine, line)
        arcAmp                  = shapely.hausdorff_distance(infLine, arcLine) # Calculate Arc Amplitde
        bendWidth, bendMaxWidth = get_bend_width(line, arcLine, dfNodeR, dfR) # get arc bendWidth and max bendWidth
        if (arcAmp / bendMaxWidth) < 0.5:
            tp.pop(1)



    
    infCoords = coords[tp] # get list of coordinates for inflection points

    lenApex = len(infCoords)-1 # Number of bends
    # create empty lists and arrays to be filled with bend values
    bendLines, infLines   = np.empty(lenApex, dtype = object), np.empty(lenApex, dtype = object)
    apexPList, apexPOList = np.empty(lenApex, dtype = object), np.empty(lenApex, dtype = object)
    bendSin   , bendWidths, bendMaxWidths       = np.empty(lenApex), np.empty(lenApex), np.empty(lenApex) 
    amplitudes, curveList , bendDO, bendLengths = np.empty(lenApex), np.empty(lenApex), np.empty(lenApex), np.empty(lenApex)
    
    ##############################
    # loop over bends/apex values
    ##############################
    for ip in range(lenApex):
        p1 = infCoords[ip]
        p2 = infCoords[ip+1]

        # get indices for curvature range
        curveIndex1 = tp[ip] if ip == 0 else tp[ip] + 1
        curveIndex2 = tp[ip+1] +1 if (ip +2) < len(tp)  else tp[ip+1]
        
        # get inflection line and bend values
        infLine = LineString([p1, p2])
        bendLine                = create_bend_line(infLine, line)
        bendWidth, bendMaxWidth = get_bend_width(line, bendLine,
                                                 dfNodeR, dfR)
    
        ##########
        # Determine apex points
        ##########
        # curvature based
        # apexP, apexPO, amplitude = get_apex_curve(bendLine, infLine, arcSign[ip])

        # distance based
        apexP, apexPO, amplitude = get_apex_distance(bendLine, infLine, arcSign[ip])
        # dt52 = dt.now()

        #####################
        # save values
        #####################
        #Geometries
        bendLines[ip]  = bendLine
        infLines[ip]   = infLine
        apexPList[ip]  = apexP
        apexPOList[ip] = apexPO

        # valies
        bendLengths[ip]   = bendLine.length
        bendSin[ip]       = bendLine.length / infLine.length 
        bendWidths[ip]    = bendWidth
        bendMaxWidths[ip] = bendMaxWidth
        amplitudes[ip]    = amplitude

        curveList[ip]     = np.nanmean(curve[curveIndex1:curveIndex2])
        bendDO[ip]        = get_bend_dist_out(line, bendLine, dfNodeR)

    # Reach sinuosity
    sin  = line.length / LineString(infCoords).length

    return (sin, bendSin,infCoords, coords[curveChanges], amplitudes, apexPList, apexPOList, curveList, 
                bendLines, infLines, bendWidths, bendMaxWidths, bendDO, bendLengths)