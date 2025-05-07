#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions:
- create_dir
- get_ids
- point_csv_to_shapely_point
"""
# import packages
import os
import shapely
import pandas as pd
import geopandas as gpd
import numpy as np
import ast
import json
import psutil
from osgeo import gdal

# import partial packages
from shapely.geometry import LineString
from glob import glob
from matplotlib.colors import LinearSegmentedColormap, ListedColormap


# import custom modules
from calc_functions import get_entrenchment_slope_intersect,\
      get_entrenchment_slope_no_intersect, slope_curvature, confinement_values
# from confinement_margin import confinement_margin


def find_connected_side(node1, node2):
    """
    Find the connected nodes to different reach. The values can be either [0,1] or [-1, -2]
    input:
    - node1: nodes of target reach
    - node2: nodes of connected reach\n
    output:
    - list with row numbers for the connection.
    """
    nodeS = node2.geometry.distance(node1.geometry.iloc[0]).min()
    nodeE = node2.geometry.distance(node1.geometry.iloc[-1]).min()
    if nodeS > nodeE:
        rowConnect = [-1, -2]
    else:
        rowConnect = [0,1]
    return rowConnect

def vw_one_side(hp, d, threshold):
    dh = list(zip(d, hp))
    dt = list(zip(d, [threshold + hp[0]] * len(d)))
    dhLine = LineString(dh)
    dtLine = LineString(dt)

    if dhLine.intersects(dtLine):
        inter = dhLine.intersection(dtLine)
        x     = shapely.get_coordinates(inter)[:,0]
        vw    = abs(x).min()
    else:
        vw = np.nan
    return vw, dhLine, dtLine

# def alternative_entrenchment_ratio(hp1List, hp2List, dInnList,width, thresholdRatio, plot, plotRow):
    """
    Find valley width and compute valley entrenchment
    input:
        - hp1List/hp2List: elevation profile too both sides of the river centerline list of lists
        - dInnList: distance list: only one needed since distance is equal to both sides list of lists
        - width: list of bendwidths
        - thresholdRatio: multiplier for the width to determine the threshold height
        - plot: boolean for plotting results
        - plotRow: dataframe row values for geometry and apex locations
    """
    
    
    
    if (isinstance(hp1List, list) == False) | (isinstance(hp1List, list) == False):

        return np.nan, np.nan, np.nan,np.nan, np.nan, np.nan, np.nan,np.nan, np.nan,np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    
    vwList, singleSideList, erS1List, erS2List, SInnList, SOutList, slopeCurveOutList, slopeCurveInnList = [], [], [], [], [], [], [], []
    erLeft, erRight, SLeft, SRight, slopeCurveLeft, slopeCurveRight = [], [], [], [], [], []
    for i in range(len(hp1List)):
        hp1  = hp1List[i]
        hp2  = hp2List[i]
        dInn = dInnList[i]

        if (all(np.isnan(hp1))) | (all(np.isnan(hp2))):
            # if all inside or all outside bend values are not present discard reach ID
            vw = np.nan
            vwList.append(vw)
            singleSideList.append(vw)
            erS1List.append(vw)
            erS2List.append(vw)
            SInnList.append(vw)
            SOutList.append(vw)
            slopeCurveOutList.append(vw)
            slopeCurveInnList.append(vw)

        else:
            pr   = np.array(hp1.copy())
            darr = np.array(dInn.copy())
            pr   = np.flip(pr, 0)
            pr   = np.append(pr, hp2)

            d = darr * -1
            d = np.flip(d, 0)
            d = np.append(d, dInn)

            threshold = width[i] * thresholdRatio # threshold = width depth relation the depth is 1/28 width

            halfWidth = width[i] / 2

            riverarea = np.where((d > (-1*halfWidth)) & (d < halfWidth))[0]
            
            prRiver = pr.copy()
            prRiver = prRiver[riverarea]
            prThreshold = prRiver.min()
            
            prminId = np.where(prRiver == prThreshold)[0]
            rivMinDist = d[riverarea[prminId]]

            A = prminId[int(len(prminId) / 2)]

            centerInt  = riverarea[A]
            centerDist = d[centerInt]



            centerDist = rivMinDist[np.abs(rivMinDist).argmin()]
            centerInt = np.where(d == centerDist)[0][0]


            dh = list(zip(d, pr))

            thresholdheight = threshold + prThreshold
            dt = list(zip(d, [thresholdheight] * len(d)))
            

            dhLine = LineString(dh)
            dtLine = LineString(dt)

            ########################################################
            # Check if valley edge intersects with height profile
            ########################################################
            if dhLine.intersection(dtLine): # Intersection with height profile True
                # check all intersections
                inter = dhLine.intersection(dtLine)
                x = shapely.get_coordinates(inter)[:,0] # get x coordinates of intersection

                # Conditional statement to check if intersection if is inside (below center) or outside bend (above center)
                cond = x > centerDist 

                ###############
                # Intersection of height profile single sided or double sided
                ###############
                if (all(cond)) | (not any(cond)): # single sided intersection

                    x1 = abs(x).min() - abs(centerDist) # distance of lowest point in river to intersection
                    x2 = np.max(dInn) - abs(centerDist) # distance of lowest point in river to end of profile

                    if np.min(x) > centerDist: # Intersection is Outside bend
                        outInter, innInter = True, False # boolean for which side intersects

                        xSlope = np.min(x) # distance to intersection
                        if xSlope < halfWidth: # check if intersection lies in inside river boundary
                            SOut          = 90
                            slopeCurveOut = 1
                            erSide1       = halfWidth
                        else:

                            SOut          = get_entrenchment_slope_intersect(d, prThreshold,thresholdheight, riverarea[-1], xSlope, 'positive')
                            slopeCurveOut = slope_curvature(d, pr, thresholdheight, centerInt, xSlope)
                            erSide1       = xSlope
                        
                        # calc ER slope/ER shape/ER distance for not intersecting profile

                        SInn = get_entrenchment_slope_no_intersect(d, pr, centerInt , np.min(d), pr[0] ,'negative')
                        slopeCurveInn = slope_curvature(d, pr, thresholdheight, centerInt, np.min(d))
                        erSide2 = x2
                        
                        outX, innX = xSlope,0
                        
                    else: # intersection if inside bend
                        outInter, innInter = False, True # boolean for which side intersects

                        xSlope = np.max(x) # distance to intersection
                        if abs(xSlope) < halfWidth: # check if intersection lies in inside river boundary
                            SInn = 90
                            erSide2 = halfWidth
                            slopeCurveInn = 1
                        else:

                            SInn = get_entrenchment_slope_intersect(d, prThreshold,thresholdheight, riverarea[0],  xSlope, 'negative')
                            slopeCurveInn = slope_curvature(d, pr, thresholdheight, centerInt, xSlope)
                            erSide2 = abs(xSlope)

                        # calc ER slope/ER shape/ER distance for not intersecting profile

                        SOut = get_entrenchment_slope_no_intersect(d, pr, centerInt, np.max(d), pr[-1],'positive')
                        slopeCurveOut = slope_curvature(d, pr, thresholdheight, centerInt, np.max(d))
                        erSide1 = x2

                        outX, innX = 0,xSlope
        
                    # get Valley width 
                    vw = x1 + x2 
                    # set boolean for single or double sided Valley confinement                   
                    confSingle = True
                    
                else: # Double sided intersection

                    x1, x2 = [x[cond], x[~cond]]
                    x1 = x1.min() # positive side
                    x2 = x2.max() # negative side
                    outInter, innInter = True, True
                    outX, innX = x1,x2

                    if abs(x1) < halfWidth:
                        x1 = halfWidth
                        SOut = 90
                        slopeCurveOut = 1
                    else: # positive side
                        # SOut = 1
                        SOut = get_entrenchment_slope_intersect(d, prThreshold,thresholdheight, riverarea[-1], x1, 'positive')

                        slopeCurveOut = slope_curvature(d, pr, thresholdheight, centerInt, x1)


                    if abs(x2) < halfWidth:
                        x2 = -1 * halfWidth
                        SInn = 90
                        slopeCurveInn = 1
                    else:

                        SInn = get_entrenchment_slope_intersect(d, prThreshold,thresholdheight, riverarea[0], x2, 'negative')
                        slopeCurveInn = slope_curvature(d, pr, thresholdheight, centerInt, x2)
                        


                    erSide1 = abs(x1)
                    erSide2 = abs(x2)
                    vw = abs(np.diff([x1,x2])[0])
                    confSingle = False
            else: # Intersection with height profile False
 
                vw = np.max(dInn) * 2
                outInter, innInter = False, False
                outX, innX = 0,0
                SOut = get_entrenchment_slope_no_intersect(d, pr, centerInt, np.max(d), pr[-1],'positive')
                SInn = get_entrenchment_slope_no_intersect(d, pr, centerInt , np.min(d), pr[0] ,'negative')
                slopeCurveOut = slope_curvature(d, pr, thresholdheight, centerInt, np.max(d))
                slopeCurveInn = slope_curvature(d, pr, thresholdheight, centerInt, np.min(d))


                erSide1, erSide2 = np.nan, np.nan
                confSingle = np.nan

                # get slope in two direction for total area
            
            
            ###############
            # Add all relevant bend values to list for entire reach
            ###############
            vwList.append(vw) # valley width
            singleSideList.append(confSingle) # single/double/no intersection
            
            # single side valley width
            erS1List.append(erSide1) # Outside bend
            erS2List.append(erSide2) # Inside bend

            # slope inside and outside bend
            SOutList.append(SOut)
            SInnList.append(SInn)

            # concavity/convexity of entrenchment
            slopeCurveOutList.append(slopeCurveOut)
            slopeCurveInnList.append(slopeCurveInn)

            line           = plotRow.geometry
            lineOffset     = line.offset_curve(2)
            inflectionLine = LineString(plotRow.infP5_0)
            apexP          = plotRow.apexP5_0


            for a in apexP:
                ap = line.interpolate(line.project(a))
                ip = shapely.ops.nearest_points(ap,inflectionLine)[1]
                apip = LineString([ap, ip])
                if apip.intersects(lineOffset):
                    erLeft          = erSide2
                    SLeft           = SInn
                    slopeCurveLeft  = slopeCurveInn
                    erRight         = erSide1
                    SRight          = SOut
                    slopeCurveRight = slopeCurveOut
                else:
                    erLeft          = erSide1
                    SLeft           = SOut
                    slopeCurveLeft  = slopeCurveOut
                    erRight         = erSide2
                    SRight          = SInn
                    slopeCurveRight = slopeCurveInn

        if (plot == True) & (~np.isnan(vw)):
            inflectionLine = LineString(plotRow.infP5_0)
            f, ax = plt.subplots(nrows =1, ncols = 2, figsize = [15,5], width_ratios=[3, 1.2])

            ax[0].plot(*dhLine.xy, color = 'blue', label = 'Height profile')
            ax[0].plot(*dtLine.xy, color = 'limegreen', linestyle = ':', label = 'Valley Threshold')

            if innInter == True:
                X0 = d[riverarea[0]-1]
                X1 = innX
                if SInn == 90:
                    X1 = X0
                ax[0].plot([X0,X1], [prThreshold, thresholdheight],  'gold', label = 'Slope Inn', linewidth = 2, zorder = 100)
            
            else:
                dInnPlot = d[0:riverarea[0]]
                
                YPlot = ((-1*np.deg2rad(SInn)) *dInnPlot)  + prThreshold
                ax[0].plot(dInnPlot, YPlot, color = 'gold')

            if outInter == True:
                X0 = d[riverarea[-1]]
                X1 = outX
                if SOut == 90:
                    X1 = X0
                ax[0].plot([X0,X1], [prThreshold, thresholdheight], 'orange', label = 'Slope Out', linewidth = 2, zorder = 100)
            else:
                dOutPlot = d[riverarea[-1]+1::]
                YPlot = np.deg2rad(SOut) * np.arange(len(dOutPlot)) + prThreshold
                ax[0].plot(dOutPlot, YPlot, color = 'orange')

            if np.isnan(vw):
                ax[0].set_title('No Confinement')
            else:
                ax[0].set_title(f'{np.round(int(vw) / int(width[i]),2)} ({int(vw)}, {int(width[i])})\
\nERInn: {np.round((abs(erSide2)) / halfWidth,2)}, EROut: {np.round(erSide1 / halfWidth,2)}\
\nSInn: {np.round(SInn,4)}, SOut: {np.round(SOut,4)}\
\nSCInn: {np.round(slopeCurveInn,3)}, SCOut: {np.round(slopeCurveOut,3)}')
            
            ax[0].vlines(0,np.min(pr), np.max(pr), linestyle = '--', color = 'lightsteelblue', label = 'Centerline')     # centerline
            ax[0].vlines(d[riverarea[0]],np.min(pr), np.max(pr), linestyle = '--', color = 'red',label = 'River border')  # riverborder
            ax[0].vlines(d[riverarea[-1]],np.min(pr), np.max(pr), linestyle = '--', color = 'red') #riverborder
            ax[0].scatter(centerDist, prThreshold , marker = '*', label = 'New center', c = 'black', s = 70, zorder = 500)
            
            ax[0].legend()
            ax[0].grid()

            ax[1].plot(*plotRow.geometry.xy, color = 'blue', label = 'River')
            ax[1].plot(*inflectionLine.xy, zorder = 50, color = 'green', label = 'Inflection line')
            apexPoints = plotRow.apexP5_0
            c = 0
            for p in range(len(apexPoints)):
                if p != i:
                    if c == 0:
                        ax[1].scatter(*apexPoints[p].xy, color = 'black', zorder = 100, s = 10, label = 'Apex')
                    else:
                        ax[1].scatter(*apexPoints[p].xy, color = 'black', zorder = 100, s = 10)
                    c +=1
                else:
                    ax[1].scatter(*apexPoints[p].xy, color = 'red', marker = '*', zorder = 100, s = 10, label = 'Shown Apex')
            
            ax[1].axis('off')
            ax[1].axis('equal')
            ax[1].legend()
            plt.show()

            
            
    

    # calculate ratios and averages for single reach values
    vw = np.array(vwList)
    entrenchment = vw / width
    entrenchmentRatio    = np.nanmean(entrenchment)
    entrenchmentRatioMax = np.nanmax(entrenchment)
    entrenchmentRatioMin = np.nanmin(entrenchment)
    entrenchmentRatioStd = np.nanstd(entrenchment)
    valleyWidth    = np.nanmean(vw)
    valleyWidthMax = np.nanmax(vw )
    valleyWidthMin = np.nanmin(vw )
    valleyWidthStd = np.nanstd(vw )

    return (entrenchmentRatio, entrenchmentRatioMax, entrenchmentRatioMin, entrenchmentRatioStd,
            valleyWidth, valleyWidthMax, valleyWidthMin, valleyWidthStd, 
            entrenchment, vw, erS1List, erS2List,
            singleSideList,
            SOutList, SInnList,
            slopeCurveOutList, slopeCurveInnList)

 
def create_dir(path):
    check = os.path.exists(path)
    if check == False:
        os.makedirs(path)
        
def get_ids(row:'gpd.GeoDataSeries'):
    """
    Insert: geopandas dataSeries row\n
    return: list of ids
    """
    
    def ids(rch_ids):
        rch_ids = str(rch_ids).strip()
        
        if (rch_ids != 'nan') & (rch_ids != 'None') & (len(rch_ids) > 0):
            id_list = list(map(float, rch_ids.split(' ')))
            id_list = [int(A) for A in id_list]
        else:
            id_list = []
        return id_list
    rch_ids = []
    rch_ids_up = ids(row.rch_ids_up)
    rch_ids.extend(rch_ids_up)
    
    rch_ids.append(row.reach_id)
    
    rch_ids_dn = ids(row.rch_ids_dn)
    rch_ids.extend(rch_ids_dn)

    return rch_ids

# def point_csv_to_shapely_point(df):
    pointColumns = []
    for c in df.columns:
        if (c[0:5] == 'apexP') | (c[0:4] == 'infP'):   
            pointColumns.append(c)
    fail = []
    for index, row in df.iterrows():
        for col in pointColumns:
            try:
                if isinstance(row[col], float) == False:  
                    point_strings = row[col].split(',')
                    points = []
                    for i in range(len(point_strings)):
                        for pat in ['[', ']', '<', '>']:
                            point_strings[i] = point_strings[i].replace(pat, '')
                        point = shapely.wkt.loads(point_strings[i])
                        points.append(point)
                else:
                    points = np.nan
                df.at[index, col] = points
            except:

                df.at[index, col] = []
    return df

# def open_mg_file(path,geometryColumns, projection):
    df = pd.read_csv(path)

    
    # decode geometry columns as strings back into shapely objects
    for c in geometryColumns:
        for index, row in df.iterrows():
            try:
                df.loc[index, c] = shapely.wkt.loads(row[c])
            except:
                df.loc[index, c] = shapely.LineString()

    # finally reconstruct geodataframe
    df = point_csv_to_shapely_point(df)
    
    gdf = gpd.GeoDataFrame(df, crs = projection)
    return gdf

# def change_str_list(df, cols, nestedLists,totalMean, totalMeanName, ERthreshold):
    def list_to_str(cell:'str'):
        cell = cell.replace('nan', '99999') 
        listVal = ast.literal_eval(cell)
        # check if slope val is single value --> not a list
        
        
        return listVal
        # check if slope val is single value --> not a list          


    ########################################################################
    # 
    ########################################################################
    loopCounter = 0
    for index, row in df.iterrows():

        for col in cols:
                # try:
                    if isinstance(row[col] ,str): # check if val is not nan

                        R = row[col]
                        vals = list_to_str(R)
                        
                        
                        if col not in nestedLists:
                            if isinstance(vals, int):
                                listVal = [vals]
                            vals = np.array(vals)
                            vals = vals.astype(float)
                            vals[vals == 99999] = np.nan

                            df.loc[index, f'm{col}']   = np.nanmean(vals)
                            df.loc[index, f'max{col}'] = np.nanmax(vals)
                        else:
                            vals = [[X] if isinstance(X, int) else X for X in vals] # change single value to list value      
                            vals = [[np.nan if val == 99999 else val for val in inner_list] for inner_list in vals] # change missing value in nested list to nan
                        df.at[index, col] = list(vals)
                    else:
                        vals = np.nan
                        if col not in nestedLists:
                            df.loc[index, [f'{col}', f'm{col}', f'max{col}']]    = np.nan
                        else:
                            df.loc[index, f'{col}']    = np.nan
                # except:

                #     df.loc[index, [f'{col}', f'm{col}', f'max{col}']]    = np.nan
        for i, nl in enumerate(totalMean):
            newR = df.loc[index, nl].values

            if (isinstance(newR[0], float))  & (isinstance(newR[1], float)):
                newVal = np.nan
            elif (isinstance(newR[0], list)) & (isinstance(newR[1], float)):
                newVal = np.nanmean(newR[0])
            elif (isinstance(newR[0], float)) & (isinstance(newR[1], list)):
                newVal = np.nanmean(newR[1])
            else:
                newVal = np.nanmean(np.concatenate(newR))


            df.loc[index, totalMeanName[i]] = newVal

        opr       = df.loc[index, 'outPR']
        ipr       = df.loc[index, 'innPR']
        d         = df.loc[index, 'CDists']
        width     = df.loc[index, 'bendWidths']
        
        dwRelation = (1/28)
        threshold = ERthreshold *dwRelation

        # Nodig voor plot
        # geometry
        # apecices 
        # correct Apex --> in ER function
        # DEM??

        # if row.reach_id == 11444300101:
        #     PP = False
        # else:
        #     PP = False

        ###############
        # ER code
        ##############
        PP = False
        ERPlotRow = df.loc[index, ['geometry','infP5_0', 'apexP5_0']]
        er2, erMax2, erMin2, erSTD2,\
            val2, valMax2, valMin2, valSTD2,\
            ERTot2, valTot2, EROut2, ERInn2,\
            ERSingleSided,\
            SOutList,SInnList,\
            SCOutList, SCInnList = alternative_entrenchment_ratio(ipr, opr, d,width, threshold, plot = PP, plotRow = ERPlotRow)

        df.loc[index, 'ER']     = er2
        df.loc[index, 'ERMax']  = erMax2
        df.loc[index, 'ERMin']  = erMin2
        df.loc[index, 'ERStd']  = erSTD2
        df.loc[index, 'VAL']    = val2
        df.loc[index, 'VALMax'] = valMax2
        df.loc[index, 'VALMin'] = valMin2
        df.loc[index, 'VALStd'] = valSTD2
        if isinstance(SInnList, float):
            SInnList = [SInnList]
            SCInnList = [SCInnList]
        if isinstance(SOutList, float):
            SOutList = [SOutList]
            SCOutList = [SCOutList]
        
        InnOutList = SInnList.copy()
        InnOutList.extend(SCOutList)
        df.loc[index, 'IOAngle'] = np.nanmean(InnOutList)
        df.loc[index, 'InnAngle'] = np.nanmean(SInnList)
        df.loc[index, 'OutAngle'] = np.nanmean(SOutList)
        
        SCInnOutList = SCInnList.copy()
        SCInnOutList.extend(SCOutList)
        df.loc[index, 'SCInnOut'] = np.nanmean(SCInnOutList)
        df.loc[index, 'SCInn'] = np.nanmean(SCInnList)
        df.loc[index, 'SCOut'] = np.nanmean(SCOutList)

        reach = df.loc[index]
        dTotal, d1, d2 = confinement_margin(reach, df_node, directory,projection = 'EPSG:3857',distanceFactor = 2, heightFactor = ERthreshold, stepsize = 50,plot =  False)


        elevation = np.nan
        if isinstance(opr, list):
            h = [o[0] for o in opr]
            elevation = np.nanmean(h)

        df.loc[index, 'elevation'] = elevation
        
        if isinstance(ERTot2, float):
            ERTot2 = [ERTot2]


        dfER1 = pd.DataFrame({'reach_id':row.reach_id, 'ER':ERTot2, 
                            'VW':valTot2,'BW':width,
                            'ERSide1':EROut2, 'ERSide2':ERInn2, 'ERSS': ERSingleSided,
                            'SInn':SInnList, 'SOut':SOutList,
                            'SCOut':SCOutList, 'SCInn':SCInnList})
        gdfER1 = gpd.GeoDataFrame(dfER1)
        if (gdfER1.shape[0] == 1) & (np.isnan(gdfER1.iloc[0].SInn) == True):
            gdfER1['geometry'] = np.nan
        else:
            gdfER1['geometry'] = row.apexP5_0
        gdfER1 = gdfER1.set_geometry('geometry')


        if loopCounter == 0:
            dfER  = dfER1
            gdfER = gdfER1
        else:
            dfER  = pd.concat([dfER, dfER1], ignore_index = True)
            gdfER = pd.concat([gdfER, gdfER1], ignore_index = True)

            
        loopCounter +=1

    return df, dfER, gdfER

#%%

# def open_all_files(save_individual, FA, F_slice,geometryCols, strListCols,nestedLists, totalMean, totalMeanName , directory, projection, GLIM, cross_dist):
    """
    open list of dataframes and concatenate files
    input:
        save_individual: save individual files
        FA: file list
        F_slice: indices for file names
        directory: file directory
    """
    for i in range(len(FA)):
        fa = FA[i][F_slice[0]:F_slice[1]]
        
        df = open_mg_file(FA[i],geometryCols, 'EPSG:3857')

        ERThreshold = 2
        df, dfER = change_str_list(df,strListCols,nestedLists, totalMean, totalMeanName, ERThreshold)

    
        df.loc[:, 'File']   = fa
        dfER.loc[:, 'File'] = fa
        df['ESSlopeRatio'] = (df['moutSlope'] - df['lineSlope']) / (((df['moutSlope'] + df['lineSlope'])/2))
        df['ERSlopeRatio'] = (df['minnSlope'] - df['lineSlope']) / (((df['minnSlope'] + df['lineSlope'])/2))
        
        ##### ADD GLIM MATCH
        
        dfG = gpd.sjoin(df, GLIM, how='left', op='intersects')
        dfG = dfG.reset_index(drop = True)

        # Calculate the geometry of the intersection
        dfG['intersection'] = dfG.apply(
            lambda row: row['geometry'].intersection(row['GLGeometry']), axis=1
            )
        dfG['interSize'] = dfG['intersection'].length
        idx = dfG.groupby('reach_id')['interSize'].idxmax()
        
        idxNoNan = idx[~idx.isna()]
        idxNan   = idx[idx.isna()]

        dfGNoL = dfG[dfG.reach_id.isin(idxNan.index)]
        dfGL   = dfG.loc[idxNoNan]
        dfG    = pd.concat([dfGNoL, dfGL])
 
        df  = dfG.drop(['IDENTITY_', 'index_right', 'Shape_Length',	'Shape_Area',	'GLGeometry',	'intersection'	,'interSize'], axis = 1)

        if i > 0:
            dfT   = pd.concat([dfT, df], axis = 0, join = 'outer', ignore_index=True)
            dfERT = pd.concat([dfERT, dfER], ignore_index=True)
        else:
            dfT   = df
            dfERT = dfER
    
        if save_individual== True:
            #add save code
            savePath = directory + f'results/single_values/{fa}_{cross_dist}_{ERThreshold}_single_values.shp'
            saveCols = ['reach_id','File', 'reach_len','width','node_mwm','n_chan_max', 'n_chan_mod','lineSlope','facc', 'dist_out','slope', 'cycle','geometry',
                        'river_name',
                        'LCheck', 'min_len', 'short',
                         'widthRatio', 
                         'moutSlope', 'minnSlope',
                         'mang5_0', 'sin5_0',
                         'ER', 'ERMax', 'ERMin', 'ERStd',
                         'VAL', 'VALMax', 'VALMin', 'VALStd',
                         'xx','Litho',
                         'elevation',
                         'IOAngle', 'InnAngle', 'OutAngle',
                         'SCInnOut']    
            saveCols.extend(totalMeanName)
            df = df[saveCols].to_file(savePath, driver='ESRI Shapefile', crs = projection)
            savePathER = directory + f'results/single_values/{fa}_{cross_dist}_ER_{ERThreshold}_values.csv'
            dfER.to_csv(savePathER)

    return dfT
#%%
# def open_all_files_multiproces(save_individual, FA, F_slice,geometryCols, strListCols,nestedLists, 
#                              totalMean, totalMeanName , directory, projection, GLIM, cross_dist, ERThreshold):
    """
    open list of dataframes and concatenate files
    input:
        save_individual: save individual files
        FA: file list
        F_slice: indices for file names
        directory: file directory
    """
    
    fa = FA[F_slice[0]:F_slice[1]]

    df = open_mg_file(FA,geometryCols, 'EPSG:3857')
    # df_node = gpd.read_file(directory + '')
    df, dfER, gdfER = change_str_list(df,strListCols,nestedLists, totalMean, totalMeanName, ERThreshold)


    df.loc[:, 'File']    = fa
    dfER.loc[:, 'File']  = fa
    gdfER.loc[:, 'File'] = fa

    df['ESSlopeRatio'] = (df['moutSlope'] - df['lineSlope']) / (((df['moutSlope'] + df['lineSlope'])/2))
    df['ERSlopeRatio'] = (df['minnSlope'] - df['lineSlope']) / (((df['minnSlope'] + df['lineSlope'])/2))
    
    ##### ADD GLIM MATCH
    
    dfG = gpd.sjoin(df, GLIM, how='left', op='intersects')
    dfG = dfG.reset_index(drop = True)

    # Calculate the geometry of the intersection
    dfG['intersection'] = dfG.apply(
        lambda row: row['geometry'].intersection(row['GLGeometry']), axis=1
        )
    dfG['interSize'] = dfG['intersection'].length
    idx = dfG.groupby('reach_id')['interSize'].idxmax()
    
    idxNoNan = idx[~idx.isna()]
    idxNan   = idx[idx.isna()]

    dfGNoL = dfG[dfG.reach_id.isin(idxNan.index)]
    dfGL   = dfG.loc[idxNoNan]
    dfG    = pd.concat([dfGNoL, dfGL])

    df  = dfG.drop(['IDENTITY_', 'index_right', 'Shape_Length',	'Shape_Area',	'GLGeometry',	'intersection'	,'interSize'], axis = 1)



    if save_individual== True:
        #add save code
        savePath = directory + f'results/single_values/{fa}_{cross_dist}_{ERThreshold}_single_values.shp'
        saveCols = ['reach_id','File', 'reach_len','width','node_mwm','n_chan_max', 'n_chan_mod','lineSlope','facc', 'dist_out','slope', 'cycle','geometry',
                    'river_name',
                    'LCheck', 'min_len', 'short',
                        'widthRatio', 
                        'moutSlope', 'minnSlope',
                        'mang5_0', 'sin5_0',
                        'ER', 'ERMax', 'ERMin', 'ERStd',
                        'VAL', 'VALMax', 'VALMin', 'VALStd',
                        'xx','Litho',
                        'elevation',
                        'IOAngle', 'InnAngle', 'OutAngle',
                        'SCInnOut']    
        saveCols.extend(totalMeanName)
        df[saveCols].to_file(savePath, driver='ESRI Shapefile', crs = projection)
        
        savePathER = directory + f'results/single_values/{fa}_{cross_dist}_ER_{ERThreshold}_values.csv'
        dfER.to_csv(savePathER)
        savePathER = directory + f'results/single_values/{fa}_{cross_dist}_ER_{ERThreshold}_values.shp'
        gdfER.to_file(savePathER, driver='ESRI Shapefile', crs = projection)
        

    return df, dfER, gdfER



def get_local_utm_projection(DF):
    """
    Calculates the local UTM projection for each row in dataframe based on x and y column in dataframe
    input:
        DF: geopandas dataframe
    output:
        Updated geopandas dataframe
    """
    df = DF.copy() # create copy to avoid linking 
    
    for I, R in df.iterrows():
        lat = R.y
        lon = R.x
        
        # Calculate the UTM zone number
        zone_number = int((lon + 180) / 6) + 1

        # 326 is north and 327 is south     
        epsg_code = 'EPSG:326' if lat >= 0 else 'EPSG:327'

        # if zone_number below 10 add 0 before value
        if zone_number < 10:
            epsg_code += f'0{zone_number}'
        else:
            epsg_code += f'{zone_number}'
        
        # add value to row in dataframe
        df.loc[I, 'localCRS'] = epsg_code
    return df

def change_crs_of_list(valList, crsCurrent, crsNew):
    """
    change crs of a list of geometries
    input:
        valList: the to change values
        crsCurrent: crs of current values
        crsNew: new crs
    output:
        new list with changed geometries
    """
    dfChange = gpd.GeoDataFrame({'geometry': valList}, crs = crsCurrent)
    dfChange = dfChange.to_crs(crsNew)
    newList = dfChange.geometry.tolist()
    return newList

def create_custom_cmap(colors : 'list',continuous = False, demMap  : 'bool'= True):
    """Create custom colormap \n
        input: \n
            colors: colors to be used in colormap. First color is lower value \n
            demMap: Select predetermined dem colormap (default = True) \n
        output: cmap
    """
    if demMap == True:
        colors = ['purple', 'blue', 'cyan', 'green', 'yellow', 'red']
    if continuous == True:
        cmap = LinearSegmentedColormap.from_list("RedYellowGreen", colors)
    else:
        cmap = ListedColormap(colors)
    return cmap

def SWORD_stats(df, directory):
    dfT = df.copy()

    dfInc = dfT[dfT['include_flag'] == '0']

    dfStat = pd.DataFrame({'continent': 'Total', 
                            'total_length'       : int(dfT['reach_len'].sum())                 , 'total_length_filtered': int(dfInc['reach_len'].sum()),
                            'num_reaches'        : dfT.shape[0]                                , 'num_reaches_filtered' : dfInc.shape[0],
                            'min_length'         : int(dfT['reach_len'].min())                 , 'max_length'           : int(dfT['reach_len'].max()),
                            'min_length_filtered': int(dfInc['reach_len'].min())               , 'max_length_filtered'  : int(dfInc['reach_len'].max()),
                            'min_length_comb'    : int(dfInc['combined_reach_len'].min())      , 'max_length_comb'      : int(dfInc['combined_reach_len'].max()),
                            'min_width_comb'     : int(dfInc['combined_reach_width'].min())    , 'max_width_comb'       : int(dfInc['combined_reach_width'].max()),
                            'width_quantile95'   : dfInc['combined_reach_width'].quantile(0.95), 'max_width_quantile95' : dfInc['combined_reach_max_width'].quantile(0.95)}, index = [0])
        

    for c in dfT['file_cont'].unique():

        cont = dfT[dfT['file_cont'] == c].copy()
        contF = cont[cont['include_flag'] == '0']
        temp = pd.DataFrame({'continent': c, 
                            'total_length'       : int(cont['reach_len'].sum())             , 'total_length_filtered': int(contF['reach_len'].sum()),
                            'num_reaches'        : cont.shape[0]                            , 'num_reaches_filtered' : contF.shape[0],
                            'min_length'         : int(cont['reach_len'].min())             , 'max_length'           : int(cont['reach_len'].max()),
                            'min_length_filtered': int(contF['reach_len'].min())           , 'max_length_filtered'  : int(contF['reach_len'].max()),
                            'min_length_comb'    : int(contF['combined_reach_len'].min())  , 'max_length_comb'      : int(contF['combined_reach_len'].max()),
                            'min_width_comb'     : int(contF['combined_reach_width'].min()), 'max_width_comb'       : int(contF['combined_reach_width'].max()),
                            'width_quantile95'   : contF['combined_reach_width'].quantile(0.95), 'max_width_quantile95' : contF['combined_reach_max_width'].quantile(0.95)}, index = [0])

        dfStat = pd.concat([dfStat, temp],ignore_index = True)
    dfStat.to_csv(directory + 'results/SWORD_stats.csv')

def shaped_logarithmic(x,quant, FMax, FMin, shape=1.0):
    
    # remove highest values (outliers)
    x[x > x.quantile(quant)] = x.quantile(quant)
    x = x.values

    xMin, xMax = x.min(), x.max()

    normalized = (np.log(x) - np.log(xMin)) / (np.log(xMax) - np.log(xMin))  # Normalize log(x)
    return FMin + (FMax - FMin) * normalized**shape  # Apply shape transformation

def shaped_lineair(x, quantL, quantH, y2, y1):
    
    # remove highest values (outliers)
    x[x > x.quantile(quantH)] = x.quantile(quantH)
    x[x < x.quantile(quantL)] = x.quantile(quantL)

    x = x.values

    xMin, xMax = x.min(), x.max()
    a = (y2 - y1) / (xMax - xMin)
    b = y1 - (((y2 - y1)*xMin) / (xMax - xMin) )

    return (a*x +b)

def smooth_factor(df, directory):
    dfInc  = df[df['include_flag'] == '0'].copy()

    dfIncG = dfInc.groupby('combined_reach_id').first()

    # dfIncG['smoothFactor'] = shaped_logarithmic(dfIncG['combined_reach_width'].copy() ,0.95, 1, 5, 1)
    dfIncG['smoothFactor'] = shaped_lineair(dfIncG['combined_reach_width'].copy() ,0.95, 1, 5)
    
    dfSM = dfIncG[['smoothFactor', 'combined_reach_width']].reset_index().copy()
    dfSM.to_csv(directory + 'results/smoothingFactor.csv')

def file_sorting(dfT, directory):
    dfSize = dfT.loc[dfT['include_flag'] == '0'].groupby(['file'],as_index = False).size()
    dfSize = dfSize.sort_values('size', ascending = False)

    files = np.sort(glob(directory + 'results/new_segments/vector/??_??_*.gpkg'))
    fNames = [f'{file[-29:-27]}_{file[-26:-24]}' for file in files]
    dfFiles = pd.DataFrame({'file': fNames, 'filePath': files})
    
    dfFiles = dfFiles.merge(dfSize, how = 'left', on = 'file').sort_values('size', ascending = False)
    dfFiles.to_csv(directory +'results/file_sorting.csv')

def node_position(line, dfN):
    dfN = dfN.copy()
    for i, r in dfN.iterrows():
        linePos = line.project(r.geometry)
        dfN.loc[i, 'linePos'] = linePos
    return dfN

def adjust_new_segments(df):
    """Change values that are saved as str but should be lists"""
    for i, r in df.iterrows():
        for t in ['', '_orig']:
            for side in ['up', 'dn']:
                rch = r[f'rch_id_{side}{t}']
                if rch == 'nan':
                    rch = np.nan
                else:
                    rch = ast.literal_eval(rch)
                df.at[i, f'rch_id_{side}{t}'] = rch
    return df

def str_to_list(df, listCols, nestedListCols):
    for i, r in df.iterrows():
        for col in listCols:
            
            colVal = r[col]
            if isinstance(colVal, float):
                colVal = np.nan
            elif isinstance(colVal, list):
                continue
            else:
                if col in nestedListCols:
                    rc = r[col].replace('(', '').replace(')', '')
                    colVal = json.loads(rc)
                else:
                    colVal = ast.literal_eval(r[col])
            df.at[i, col] = colVal
    return df

def str_to_list_comb(df, listCols, nestedListCols):

    dfN = df.groupby('combined_reach_id', as_index = False).first()

    for i, r in dfN.iterrows():
        reaches = df.loc[df['combined_reach_id'] == r['combined_reach_id'], 'reach_id'].values
        dfN.at[i , 'reach_ids'] = str(reaches)
        for col in listCols:
            colVal = r[col]

            if isinstance(colVal, float):
                colVal = np.nan
            elif isinstance(colVal, list):
                continue
            else:
                if col in nestedListCols:
                    rc = r[col].replace('(', '').replace(')', '')
                    colVal = json.loads(rc)
                else:
                    colVal = ast.literal_eval(r[col])
            dfN.at[i, col] = colVal
    return dfN

def expand_dataframe(df):
    """
    Expands a DataFrame by creating new rows for columns containing lists,
    while repeating single-valued columns.
    
    :param df: pandas DataFrame
    :return: Expanded pandas DataFrame
    """
    list_cols = [col for col in df.columns if df[col].apply(lambda x: isinstance(x, list)).any()]

    if not list_cols:
        return df  # Return original if no list columns found
    
    # Expand DataFrame
    df_expanded = df.explode(list_cols, ignore_index=True)
    return df_expanded

def confinement_factor(df, y1, y2):

    dfIncG = df.groupby('combined_reach_id', as_index = False).first()
    cf = shaped_lineair(dfIncG['combined_reach_width'].copy() ,0.05,0.95, y1, y2)
    
    dfIncG['conFactor'] = cf = 1/cf

    dfIncG = dfIncG[['combined_reach_width', 'conFactor']]
    return dfIncG

def confinement_factor_single_values(df, widthCol, y1, y2):

    cf = shaped_lineair(df[widthCol].copy() ,0.05,0.95, y1, y2)
    
    df['conFactor'] = cf = 1/cf

    df = df[[widthCol, 'conFactor']]
    return df

def check_memory():
    process = psutil.Process(os.getpid())
    # print(f"Memory usage: {process.memory_info().rss / (1024 * 1024):.2f} MB")
    return np.round(process.memory_info().rss / (1024 * 1024), 2)

def concat_csv_files(files):
    for i, f in enumerate(files):
        dfTemp = pd.read_csv(f, dtype = {'include_flag':str, 'calculated':str})
        if i == 0:
            df = dfTemp
        else:
            df = pd.concat([df, dfTemp])
    return df

def concat_gpkg_files(files):
    for i, f in enumerate(files):
        dfTemp = gpd.read_file(f)
        if i == 0:
            df = dfTemp
        else:
            df = pd.concat([df, dfTemp])
    return df

def open_single_values(files):
    for i, f in enumerate(files):
        dfTemp = pd.read_csv(f, dtype = {'include_flag':str, 'calculated':str})
        dfTemp = str_to_list(dfTemp, ['distOut', 'distInn', 'elevInn', 'elevOut'], ['distOut', 'distInn', 'elevInn', 'elevOut'])
        if i == 0:
            df = dfTemp
        else:
            df = pd.concat([df, dfTemp])
    return df



# # Function for properly constrained exponential curves
# def shaped_exponential(x, FMin, FMax, shape=1.0):
#     x[x > x.quantile(0.95)] = x.quantile(0.95)
#     x = x.values
#     xMin, xMax = x.min(), x.max()

#     b = np.log(FMax) / (xMax - xMin)  # Ensures correct scaling
#     normalized = (np.exp(b * (x - xMin)) - 1) / (np.exp(b * (xMax - xMin)) - FMin)  # Normalize 0-1
#     return FMin + (FMax - FMin) * normalized**shape  # Scale to [1,5]

# def shaped_linear(x, FMin, FMax):
#     x[x > x.quantile(0.95)] = x.quantile(0.95)
#     x = x.values
#     m = (FMax - FMin) / (x.max() - x.min())
#     b = FMin - (m * x.min())
#     return (m * x + b)