#import packages
import numpy as np
import shapely
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import ast

#import seperate modules
from shapely import LineString

# Import custom modules
from calc_functions import get_entrenchment_slope_intersect, get_entrenchment_slope_no_intersect, slope_curvature
from confinement_margin import confinement_margin
from dem import find_dem

def point_csv_to_shapely_point(df):
    # pointColumns = []
    # for c in df.columns:
    #     # if (c[0:5] == 'apexP') | (c[0:4] == 'infP'):   
    #     if  (c[0:4] == 'infP'):   
    #         pointColumns.append(c)
    # fail = []
    # for index, row in df.iterrows():
    #     for col in pointColumns:
    #         try:
    #             if isinstance(row[col], float) == False:  
    #                 point_strings = row[col].split(',')
    #                 points = []
    #                 for i in range(len(point_strings)):
    #                     for pat in ['[', ']', '<', '>']:
    #                         point_strings[i] = point_strings[i].replace(pat, '')
    #                     point = shapely.wkt.loads(point_strings[i])
    #                     points.append(point)
    #             else:
    #                 points = np.nan
    #             df.at[index, col] = points
    #         except:

    #             df.at[index, col] = []
    
    wktCol = ['bendLines', 'apexP5_0', 'infP5_0']
    for index, row in df.iterrows():
        for col in wktCol:
            try:
                if isinstance(row[col], float) == False: 

                    colList = eval(row[col])
                    geom = [shapely.wkt.loads(A) for A in colList]
                else:
                    geom = np.nan
                df.at[index, col] = geom
            except:
                df.at[index, col] = np.nan

    return df

def alternative_entrenchment_ratio(hp1List, hp2List, dInnList,width, thresholdRatio, directory, plot, dfplotRow):

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
    plotRow = dfplotRow.iloc[0]
    
    
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

            #########################
            # right plot
            demWidth = plotRow.node_mwm * 3
            raster = find_dem(dfplotRow, directory, 'EPSG:3857',
                        demWidth*5, False) # *1.2 to create extra buffer
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


def open_mg_file(path,geometryColumns, projection):
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

def change_str_list(df, df_node, cols, nestedLists,totalMean, totalMeanName, ERthreshold, directory):
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
        # print(row.reach_id, index, loopCounter)
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

        ###############
        # ER code
        ##############
        opr       = df.loc[index, 'outPR']
        ipr       = df.loc[index, 'innPR']
        d         = df.loc[index, 'CDists']
        width     = df.loc[index, 'bendWidthM']

        dwRelation = (1/28)
        threshold = ERthreshold *dwRelation


        PP = False
        ERPlotRow = df.loc[df.index == index, ['reach_id','node_mwm', 'geometry','infP5_0', 'apexP5_0']]
        er2, erMax2, erMin2, erSTD2,\
            val2, valMax2, valMin2, valSTD2,\
            ERTot2, valTot2, EROut2, ERInn2,\
            ERSingleSided,\
            SOutList,SInnList,\
            SCOutList, SCInnList = alternative_entrenchment_ratio(ipr, opr, d,width, threshold,directory, plot = PP, dfplotRow = ERPlotRow)

        df.loc[index, 'ER']     = er2
        df.loc[index, 'ERMax']  = erMax2
        df.loc[index, 'ERMin']  = erMin2
        df.loc[index, 'ERStd']  = erSTD2
        df.loc[index, 'VAL']    = val2
        df.loc[index, 'VALMax'] = valMax2
        df.loc[index, 'VALMin'] = valMin2
        df.loc[index, 'VALStd'] = valSTD2
        
        if isinstance(SInnList, float):
            SInnListM, SCInnListM  = [SInnList], [SCInnList]
        else:
            SInnListM, SCInnListM = SInnList, SCInnList
        if isinstance(SOutList, float):
            SOutListM, SCOutListM  = [SOutList], [SCOutList]
        else:
            SOutListM, SCOutListM = SOutList, SCOutList

        InnOutList = SInnListM.copy()
        InnOutList.extend(SCOutListM)
        df.loc[index, 'IOAngle'] = np.nanmean(InnOutList)
        df.loc[index, 'InnAngle'] = np.nanmean(SInnListM)
        df.loc[index, 'OutAngle'] = np.nanmean(SOutListM)
        
        SCInnOutList = SCInnListM.copy()
        SCInnOutList.extend(SCOutListM)
        df.loc[index, 'SCInnOut'] = np.nanmean(SCInnOutList)
        df.loc[index, 'SCInn']    = np.nanmean(SCInnListM)
        df.loc[index, 'SCOut']    = np.nanmean(SCOutListM)

        reach = df.loc[df.index == index]
        if isinstance(reach.iloc[0].bendLines, list):
            dTotal, d1, d2 = confinement_margin(reach, df_node, directory,projection = 'EPSG:3857',distanceFactor = 2, heightFactor = ERthreshold, stepsize = 50,plot =  False)
        else:
            dTotal, d1, d2 = np.nan, np.nan, np.nan

        elevation = np.nan
        if isinstance(opr, list):
            h = [o[0] for o in opr]
            elevation = np.nanmean(h)

        df.loc[index, 'elevation'] = elevation
        


        ############################
        # Check for NAN values!!!!!!!!!!!!! wkt not going to work and neither are the list comprehensions 
        ############################
        infP            = df.loc[index, 'infP5_0']
        apexP           = df.loc[index, 'apexP5_0']
        ang             = df.loc[index, 'ang5_0']
        bendLines       = df.loc[index, 'bendLines']
        widthW          = df.loc[index, 'bendWidths']

        if (isinstance(infP, float)) | (isinstance(bendLines, float)):
            infLines, widthRatio, SSin, lengthPart, apexPWKT = np.nan, np.nan, np.nan, np.nan, np.nan
            dTotal, d1, d2 = np.nan, np.nan, np.nan
            apexNumbers, apexID = np.nan, f'{row.reach_id}_0'

        else:
            infLines        = [LineString([infP[i -1], infP[i]]) for i in range(1, len(infP))]
            bendLinesLength = [I.length for I in bendLines]
            widthRatio      = np.array(widthW) / np.array(width)
            SSin            = [bl / infLines[i].length for i, bl in enumerate(bendLinesLength)]
            
            # to save linesttrings and points convert to wkt 
            infLines   = [I.wkt for I in infLines]
            bendLines  = [I.wkt for I in bendLines]
            apexPWKT   = [I.wkt for I in apexP]

            totLength  = np.sum(bendLinesLength)
            lengthPart = [bl / totLength for bl in bendLinesLength]


            reach = df.loc[df.index == index]

            dTotal, d1, d2 = confinement_margin(reach, df_node, directory,projection = 'EPSG:3857',distanceFactor = 2, 
                                                heightFactor = ERthreshold, stepsize = 50,plot =  False)
            apexNumbers = np.arange(len(apexP))
            apexID      = [f'{row.reach_id}_{an}' for an in apexNumbers]

        apexNumbers
        dfVals = {'reach_id':row.reach_id,'apexNumber':apexNumbers, 'apexID': apexID,
                            'ER':ERTot2, 'VW':valTot2,
                            'ERSide1':EROut2, 'ERSide2':ERInn2, 'ERSS': ERSingleSided,
                            'SInn':SInnList, 'SOut':SOutList,
                            'SCOut':SCOutList, 'SCInn':SCInnList, 
                            'infLines': infLines, 'apexP':apexPWKT, 'ang':ang, 'bendLine':bendLines, 
                            'width':widthW,'widthM':width,'SWidthRat':widthRatio, 
                            'SSin': SSin, 'bendPerc': lengthPart, 'geometry':apexP,
                            'confPerc': dTotal, 'confPerc1':d1, 'confPerc2':d2}
        listCheck = [len(a) for a in dfVals.values() if isinstance(a, list)]

        
        if any(listCheck) == False: # No list value in dataframe dictionary
            dfVals['reach_id'] = [dfVals['reach_id']]
        elif len(set(listCheck)) > 1:
            # print('check')
            for key, value in dfVals.items():
                if key != 'reach_id':
                    dfVals[key] = np.nan
                else:
                    dfVals[key] = [dfVals[key]]
            dfVals


        dfER1 = pd.DataFrame(dfVals)

        if loopCounter == 0:
            dfER  = dfER1
        else:
            dfER  = pd.concat([dfER, dfER1], ignore_index = True)


            
        loopCounter +=1

    # merge df and dfER for sinuosity, slope, angle, braiding index
    dfM = df[['reach_id', 'facc', 'width', 'facc', 'n_chan_max', 'n_chan_mod', 'slope', 'dist_out', 'max_width', 'cycle', 'node_mwm', 'widthRatio', 'sin5_0']]
    dfM.columns = 'R_' + dfM.columns.values

    dfER = dfER.merge(dfM, how = 'left', left_on = 'reach_id', right_on = 'R_reach_id', suffixes = [None, 'Reach'])
    dfER = dfER.drop('R_reach_id', axis = 1)
    
    dfER = dfER.set_geometry('geometry')
    return df, dfER



def open_all_files_multiproces(save_individual, FA, F_slice,geometryCols, strListCols,nestedLists, 
                               totalMean, totalMeanName , directory, projection, GLIM, cross_dist, ERThreshold):
    """
    open list of dataframes and concatenate files
    input:
        save_individual: save individual files
        FA: file list
        F_slice: indices for file names
        directory: file directory
    """
    from glob import glob
    fa = FA[F_slice[0]:F_slice[1]]
    
    
    df = open_mg_file(FA,geometryCols, 'EPSG:3857')



    nodeFile = glob(directory + f'results/new_segments/node/*{fa[0:2]}*{fa[-2::]}*.shp')[0]
    df_node  = gpd.read_file(nodeFile)

    df, dfER = change_str_list(df, df_node,strListCols,nestedLists, totalMean, totalMeanName, ERThreshold, directory)


    df.loc[:, 'File']    = fa
    dfER.loc[:, 'File']  = fa
    # gdfER.loc[:, 'File'] = fa

    df['ESSlopeRatio'] = (df['moutSlope'] - df['lineSlope']) / (((df['moutSlope'] + df['lineSlope'])/2))
    df['ERSlopeRatio'] = (df['minnSlope'] - df['lineSlope']) / (((df['minnSlope'] + df['lineSlope'])/2))
    
    ##### GLIM MATCH
    dfG = gpd.sjoin(df, GLIM, how='left', op='intersects')
    dfG = dfG.reset_index(drop = True)

    ## GLIM MATCH gdfER
    dfGER = gpd.sjoin(dfER, GLIM, how='left', op='intersects')
    dfERGlim = dfGER.reset_index(drop = True)
    dfERGlim = dfERGlim.drop_duplicates(subset = ['apexID', 'Litho'], keep = 'first', ignore_index = True)
    print(dfER.shape[0], dfERGlim.shape[0])


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
        dfERGlim.to_csv(savePathER, index = False)
        # savePathER = directory + f'results/single_values/{fa}_{cross_dist}_ER_{ERThreshold}_values.shp'
        # gdfER.to_file(savePathER, driver='ESRI Shapefile', crs = projection)
        

    return df, dfER


