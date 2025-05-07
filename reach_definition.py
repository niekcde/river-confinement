##########################################
# set directory
##########################################
directory = '/scratch/6256481/'

##########################################
# import packages
##########################################
import geopandas as gpd
import pandas as pd
import numpy as np
import shapely
import networkx as nx 
import itertools 
import os

##########################################
# partial modules
##########################################
from shapely.geometry import Point
from tqdm import tqdm

##########################################
# function typecasting
##########################################
from typing import List, Tuple

##########################################
# import custom modules
##########################################
from support import get_local_utm_projection, find_connected_side
from line_functions import azimuth

##########################################
# custom functions
##########################################
def remove_updn_connection(df, dn_connected = False):
    """
    Remove up and down connection of rows that do not have include-flag == 0
    """

    df = df.copy()
    remove = df[df['include_flag'] != '0']
    for rind, r in remove.iterrows():
        for side, contra in [['up', 'dn'], ['dn', 'up']]:
            connection = r[f'rch_id_{side}']
            if isinstance(connection, list):
                dfCon = df[df['reach_id'].isin(connection)]

                for conInd, conRow in dfCon.iterrows():
                    conSideID = conRow[f"rch_id_{contra}"]
                    if isinstance(conSideID, float):
                        continue

                    # update dn reach and number
                    conSideID = conSideID.copy()
                    conSideID.remove(r.reach_id)
                    conSideN  = conRow[f"n_rch_{contra}"]
                    conSideN -= 1
                    
                    # update dn_connected_reach                    
                    if (dn_connected == True) & (contra == 'up'):
                        dnCon = conRow['dn_connected_reach']
                        if dnCon == r['reach_id']:
                            df.loc[conInd, 'dn_connected_reach'] = np.nan

                    # update rch_id_dn/rch_id_up and n_rch_dn/n_rch_dn
                    if conSideN == 0:
                        conSideID = np.nan
                    df.at[conInd, f"rch_id_{contra}"] = conSideID
                    df.at[conInd, f"n_rch_{contra}"]  = conSideN

            # remove connection values from not included rows
            df.at[rind, f"rch_id_{side}"] = np.nan
            df.at[rind, f"n_rch_{side}"]  = 0
                    
        
    return df

def update_include_flag(df, reaches, newVal):
    df = df.copy()
    df['include_flag'] = np.where(df['reach_id'].isin(reaches), 
                            np.where(df['include_flag'] == '0', 
                                        newVal, 
                                        df["include_flag"]+ newVal),
                            df['include_flag'])  # Keep the same if condition not met
    return df

def filter_SWORD_input(df:gpd.GeoDataFrame, dfNode:gpd.GeoDataFrame) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Adjust original SWORD input file:
    - extract pfafstetter level 1 and level 2 as seperate columns
    - add local UTM projection per reach
    - change rch_id_dn and rch_id_up from str to list value (missing is NaN)
    - remove type `ghost` (6)
    - remove lakeflag `Canal` (2)
    - remove lakeflag `lake` (1) reaches longer than 7500 meters. Cut-off based on histogram  of lake reach lengths
    - add node_mean_max_width value: the mean of the max_width values of nodes for a reach. if all values are missing for a reach set width as node_mean_max_width.
    input:
    - df: original SWORD df

    output: filtered reach and node dataframe
    """
    # create copies of input to avoid chaining
    df     = df.copy()
    dfNode = dfNode.copy()
    
    df['include_flag'] = '0'

    # max lakeflagged reach length --> based on histogram of lakelengths
    maxLakeLength = 7500

    ############################################
    # Change width values to nan. If no width values exist flag 2
    ############################################
    minWidthVal = 20
    def classify_column(group):
        is_all_pos = (group > minWidthVal).all()
        is_all_neg = (group < minWidthVal).all()
        
        if is_all_pos:
            return 'positive'
        elif is_all_neg:
            return 'negative'
        else:
            return 'mixed'

    # Group by 'reach_id' and apply classification separately for width and max_width
    widthGroupClassification = dfNode.groupby('reach_id', as_index=False).agg({
                                                                        'width'    : classify_column,
                                                                        'max_width': classify_column
                                                                            }).reset_index()

    # flag reaches with completely missing width and max width (missing = negative values)
    missingWidthTotal = widthGroupClassification[(widthGroupClassification['width'] == 'negative')]
    df.loc[(df['reach_id'].isin(missingWidthTotal['reach_id'])) |
           (df['width'] < minWidthVal)                                     , 'include_flag'] = '1' # flag missing rows

    excluded = df.loc[df['include_flag'] != '0', 'reach_id'].values

    # change width values for partially width values on node basis
    dfNode.loc[(~dfNode.reach_id.isin(excluded)) & 
               ((dfNode['width'] >= minWidthVal) & (dfNode['max_width'] <= minWidthVal)), 'max_width'] = dfNode['width']
    
    dfNode.loc[(~dfNode.reach_id.isin(excluded)) & 
               ((dfNode['width'] <= minWidthVal) & (dfNode['max_width'] >= minWidthVal)), 'width']     = dfNode['max_width']

    dfNode.loc[(~dfNode.reach_id.isin(excluded)) & 
               ((dfNode['width'] < minWidthVal) & (dfNode['max_width'] < minWidthVal)), ['width', 'max_width']] = np.nan
    

    # change missing max_width values to width for reach level dataframe
    df.loc[(df['include_flag'] == '0') &
           ((df['width'] >= minWidthVal) & (df['max_width'] <= minWidthVal)), 'max_width']            = df['width']

    # change max width values to width if max width lower than width
    df.loc[df['max_width'] < df['width'], 'max_width']                = df['width']
    dfNode.loc[dfNode['max_width'] < dfNode['width'], 'max_width']    = dfNode['width']

    ############################################
    # extract pfaf 1 and two from SWORD ID
    ############################################
    df['reach_id'] = df['reach_id'].astype(str)
    df['pfaf1']    = df['reach_id'].str.slice(0,1).astype(int)
    df['pfaf2']    = df['reach_id'].str.slice(1,2).astype(int)
    df['reach_id'] = df['reach_id'].astype(int)
    df = df.copy()
    
    ############################################
    # add column with local UTM projection for each reach
    ############################################
    df = get_local_utm_projection(df)

    ############################################
    # change rch_id_up and dn to list values
    ############################################
    for i, r in df.iterrows():
        for side in ['up', 'dn']:
            rch = r[f'rch_id_{side}']
            if len(rch) > 0:
                val = list(map(int, r[f'rch_id_{side}'].split()))
            else:
                val = np.nan
            df.at[i, f'rch_id_{side}'] = val
    # set original rch_id_up and dn for future checks
    df['rch_id_up_orig'] = df['rch_id_up']
    df['rch_id_dn_orig'] = df['rch_id_dn']
    df['n_rch_up_orig']  = df['n_rch_up']
    df['n_rch_dn_orig']  = df['n_rch_dn']
    
    ############################################
    # set include flag for lakeflagged reaches with certain length and remove lakeflag 2 and type 6
    ############################################
    changeReaches = df.loc[((df['lakeflag'] == 1) & (df['reach_len'] > maxLakeLength)) | (df['lakeflag']== 2) | (df['type']== 6), 'reach_id'].values
    df            = update_include_flag(df, changeReaches, '2')
    
    ############################################
    # create node mean max width
    ############################################
    NMW = dfNode.copy()
    NMW = NMW.groupby('reach_id', as_index= False)['max_width'].mean()
    NMW = NMW.rename(columns={"max_width": "node_mean_max_width"})
    df = pd.merge(df, NMW, on ='reach_id',how='left')
    df.loc[df['node_mean_max_width'].isna(), 'node_mean_max_width'] = df['width']

    ############################################
    # remove up down flagged connections and make sure same reaches in Node and df
    ############################################
    df     = remove_updn_connection(df) # remove flagged rows from up and down connections
    dfNode = dfNode[dfNode.reach_id.isin(df['reach_id'].values)]

    # return reach and node dataframe
    return df, dfNode

def get_connection_nodes(df:gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Create junction nodes with ID for the up stream and downstream side of a reach.\n
    input:
    - df: dataframe with SWORD reaches\n
    output:
    - dfNode: dataframe with junction node_id and reach id
    """
    dfNodeCheck = pd.DataFrame(columns = ['Node_id', 'reach_id', 'up', 'dn'])
    def find_nodes(side, connection,r, df):
        if side == 'up':
            contra = 'dn'
        elif side == 'dn':
            contra = 'up'
        else:
            print('find_nodes in get_connection_nodes: wrong input value. Either up or dn!')

        
        if (r['reach_id'] in df[contra].values):             
            existingNodeRow = df.index[df[contra] == r['reach_id']].tolist()
            existingNodeVal = df.loc[existingNodeRow]['Node_id'].values[0]
            upNodeID        = existingNodeVal


            for u in connection:
                IND  = df.index.max() + 1 if not df.empty else 0
                df.loc[IND, ['Node_id', 'reach_id', side]] = [upNodeID, r['reach_id'], u]
        else:
            for u in connection:
                IND       = df.index.max() + 1 if not df.empty else 0
                upNodeID  = df['Node_id'].max() + 1 if not df.empty else 0


                if   (set(connection) & set(df[side].values)):
                    upNodeID = df[df[side].isin(connection)].Node_id.values[0]
                elif r['reach_id'] in df['reach_id'].values:
                    existingRows = df[df.reach_id == r['reach_id']]
                    if all(existingRows[side].isna()) == False:
                        upNodeID = existingRows[~existingRows[side].isna()]['Node_id'].values[0]
                
                
                df.loc[IND, ['Node_id', 'reach_id', side]] = [upNodeID, r['reach_id'], u]

        return df
    

    for i, r in df.iterrows():
        up = r.rch_id_up
        dn = r.rch_id_dn

        if isinstance(up, list):
            if len(up) > 0:
                dfNodeCheck = find_nodes('up', up,r, dfNodeCheck)
        if isinstance(dn, list):
            if len(dn) > 0:
                dfNodeCheck = find_nodes('dn', dn,r, dfNodeCheck)

    dfNode = dfNodeCheck[~dfNodeCheck.duplicated(['Node_id', 'reach_id'])].copy()
    dfNode[['up', 'dn']] = dfNode[['up', 'dn']].notna()
    return dfNode

def angle_difference(a1:float, a2:float) -> float:
    """
    Calculate angle between two orientations.\n
    input:
    - a1/a2: angles in degrees (0-360)\n
    output:
    - diff: difference in degree between input angles, with max of 180 degrees
    """
    diff = abs(a1 - a2)
    if diff > 180:
        diff = 360 - diff
    return diff

def find_connected_side_short_reach(row:'pd.Series', nodes:'gpd.GeoDataFrame', targetNode:shapely.geometry.Point) -> Tuple[float, shapely.geometry.Point, shapely.geometry.Point]:
    """
    find connected side of reach when reach has only single node.\n
    input:
    - row: selected row
    - nodes: selected reach nodes
    - targetNode: target node of connected reach\n
    output:
    - nodewidth: mean width (wet) value
    - cp1: point on centerline closest to target node
    - cp2: point on centerline furthest from target node
    """
    nodeWidth = nodes.width.mean() # average width of input nodes

    # get points at start and end of centerline
    p1 = Point(row.geometry.coords[0])
    p2 = Point(row.geometry.coords[-1])

    # determine if en or start of centerline is connected to target reach
    if p1.distance(targetNode) < p2.distance(targetNode):
        cp1, cp2 = p1, p2
    else:
        cp1, cp2 = p2, p1
    return nodeWidth, cp1, cp2

def check_manual_added_width(nodes:'gpd.GeoDataFrame', side:list) -> Tuple[list,int]:
    '''
    Check if one of the connected nodes is manually added. 
    If node is mamually added, check if one of the next two nodes is not manually added. 
    If less than two are not manually added keep original input\n
    input:
    - nodes: reach nodes to check
    - side: which side of the reach to check\n
    output:
    - side: updated connected side
    - maxNeededNodes: number of nodes needed in reach
    '''
    changed = 0 

    # check if any of the two connected nodes are "manually" added and if there are at least 4 nodes in the reach
    if (any(nodes.iloc[side].manual_add.values == 1) &\
        (nodes.shape[0] >= 4)):

        nodeSide    = nodes.iloc[side].copy().reset_index()
        sideNew     = side.copy()
        
        remove_side = nodeSide.index[nodeSide['manual_add'] == 1].tolist() # get indices of manually added nodes
        for index in sorted(remove_side, reverse=True): # loop over manually added nodes 
            del sideNew[index] # delete manually node

        # check the nearest two nodes (node 3 and 4 from the connected side)
        for i in range(2):
            if len(sideNew) == 2: # if length is two stop iteration
                continue
            
            # check if potential new value is manually added
            newSideVal = side[1] + np.diff(side)[0]
            if nodes.iloc[newSideVal].manual_add == 0:
                sideNew.append(newSideVal) # append node index for not manually added node
                changed +=1
        if len(sideNew) == 2:
            side = sideNew
    
    # the max number of nodes needed in 
    maxNeededNodes = changed + 2

    return side, maxNeededNodes

def round_geometry(geom, decimals=6):
    """
    round_geometry changes the precision of the linestring without altering the topology\n
    input:
    - geom: geometry to be changed
    - decimals: precision is changed (10**-decimal)\n
    output:
    - simplified geometry
    """
    return geom.simplify(10**-decimals, preserve_topology=True)

def wa_connection(conRow:gpd.GeoDataFrame, conNodes:gpd.GeoDataFrame, row:pd.Series, nodes:gpd.GeoDataFrame) -> Tuple[float, float]:
    """
    Calculate angle between orientation of connected reaches, based on the final two nodes in both reaches. 
    With fewer than 2 nodes in a reach, the end coordinates of the reach are used.\n
    input:
    - conRow: connected row
    - conNodes: connected node row
    - row: target tow
    - nodes: target node row\n
    output:
    - da1: angle difference in degrees
    - width: width (wet) difference between the final two nodes
    """

    # find connected side of two target node dataframes
    side    = find_connected_side(nodes, conNodes)
    conSide = find_connected_side(conNodes, nodes)
    
    # both reaches contain less than two nodes
    if (nodes.shape[0] < 2) & (conNodes.shape[0] < 2):
        # use the endpoints of the reaches two find if the first or final coordinates in the centerline are connected     
        l1 = round_geometry(row.geometry)
        l2 = round_geometry(conRow.geometry)
        nearestPR, nearestPC = shapely.ops.nearest_points(l1, l2)


        tp1 = Point(row.geometry.coords[0])
        tp2 = Point(row.geometry.coords[-1])
        if tp1.distance(nearestPR) > 10:
            p1, p2 = tp2, tp1
        else:
            p1, p2 = tp1, tp2
            
        tp1 = Point(conRow.geometry.coords[0])
        tp2 = Point(conRow.geometry.coords[-1])
        if tp1.distance(nearestPC) > 10:
            cp2, cp1 = tp2, tp1
        else:
            cp1, cp2 = tp1, tp2

        nodeWidth    = nodes.width.mean()
        conNodeWidth = conNodes.width.mean()

    else:
        # find connected side for cases with short reaches and normal reaches
        if nodes.shape[0] < 2:
            nodeWidth, p1, p2 = find_connected_side_short_reach(row, nodes, conNodes.iloc[conSide[0]].geometry)

        else:
            p1, p2    = nodes.iloc[side[0]].geometry, nodes.iloc[side[1]].geometry
            side, sideChanged       = check_manual_added_width(nodes, side)
            if nodes.shape[0] < sideChanged:
                nodeWidth = nodes['width'].mean()
            else:
                nodeWidth = nodes.iloc[side]['width'].mean()
        
        if conNodes.shape[0] < 2:
            conNodeWidth, cp1, cp2 = find_connected_side_short_reach(conRow, conNodes, nodes.iloc[side[0]].geometry)
        else:
            cp1, cp2     = conNodes.iloc[conSide[0]].geometry, conNodes.iloc[conSide[1]].geometry
            
            conSide, conSideChanged = check_manual_added_width(conNodes, conSide)
            if conNodes.shape[0] < conSideChanged:
                conNodeWidth = conNodes['width'].mean()
            else:
                conNodeWidth = conNodes.iloc[conSide]['width'].mean()
    # calculate orientation angle of the reaches
    angle        = azimuth(p2, p1)
    conAngle1    = azimuth(cp1, cp2)

    da1  = angle_difference(conAngle1, angle) # difference in orientation between the two reaches

    return da1, abs(conNodeWidth - nodeWidth)

def determine_connected_reach(row:pd.Series,df:gpd.GeoDataFrame, dfNode:gpd.GeoDataFrame, junctionReaches:gpd.GeoDataFrame) -> Tuple[int|float]:
    """Determine correct connection of reaches based on width and angle between reaches.\n
    input:
    - row: selected row
    - df: complete dataframe
    - dfNode: complete node dataframe
    - junctionReaches: dataframe with all rows connected to "row"\n
    output: connected reach_id
    - conReach: connected reach
    """
    
    ID                   = row.reach_id
    
    # set primary selection criteria. Is either minimal angle difference or a width factor differnce
    min_angle_difference = 10 # at least 10 degrees difference
    min_width_difference = 2

    # Select up and dn stream reaches connected to the junction node
    T = junctionReaches[junctionReaches['up'] == True].reach_id.unique()
    V = junctionReaches[junctionReaches['dn'] == True].reach_id.unique()

    # change junction reaches to list of reach_id's
    junctionReaches = junctionReaches.reach_id.unique()
    combinations = np.array(list(itertools.product(T, V))) # combinations between up and dn stream reaches

    # select all connected rows and nodes and change crs to local reach crs
    connectedRows     = df[df.reach_id.isin(junctionReaches)]
    connectedRows     = connectedRows.to_crs(row.localCRS)
    connectedNodeRows = dfNode[dfNode.reach_id.isin(junctionReaches)]
    connectedNodeRows = connectedNodeRows.to_crs(row.localCRS)


    # set values
    angles, widths, sameName = np.zeros(len(combinations)), np.zeros(len(combinations)), []
    count = 0
    for tid in T:
        # select target row and nodes
        tr = connectedRows[connectedRows['reach_id'] == tid].iloc[0]
        tn = connectedNodeRows[connectedNodeRows['reach_id'] == tid]
        
        for vid in V:
            # select connected row and nodes
            conRow   = connectedRows[connectedRows.reach_id == vid].iloc[0]
            conNodes = connectedNodeRows[connectedNodeRows.reach_id == vid]
            
            # check if combination has the samen river name
            if (conRow.river_name ==  tr.river_name) & (conRow.river_name != 'NODATA'):
                sameName.append(True)
            else:
                sameName.append(False)

            # calculate orientation angle between the rwo reaches
            ang, w = wa_connection(conRow, conNodes, tr, tn)

            angles[count] = ang
            widths[count] = w
            count += 1
    
    # keep the reache combinations with the same name
    if any(sameName):
        angles = angles[sameName]
        widths = widths[sameName]
        combinations = combinations[sameName]
    
    # sorted angle and width indices
    sortedAngles = np.argsort(angles)
    sortedWidths = np.argsort(widths)

    # if only one combination left select this combination ELSE check width and angle difference
    if len(combinations) == 1:
        conReach = combinations
    else:
        # if abs(angles[sortedAngles[0]] - angles[sortedAngles[1]]) < min_angle_difference: # angle must be larger than threshold
        #     if idMinAngle != idMinWidth: # if width corresponds to angle keep angle selected connection otherwise follow width
        #         conReach = combinations[sortedWidths[0]]
        #     else:
        #         conReach = combinations[sortedAngles[0]]
        
        # if the smallest width difference is zero keep that reach combination
        if (widths[sortedWidths[0]] == 0) & (widths[sortedWidths[1]] != 0):
            conReach = combinations[sortedWidths[0]]
        # if both width differences are non zero
        elif (widths[sortedWidths[0]] != 0) & (widths[sortedWidths[1]] != 0):
            # check if ratio between the numbers is below minimal preset ratio difference
            if (widths[sortedWidths[1]] / widths[sortedWidths[0]]) < min_width_difference:
                # if ratio is below minimal value keep angle value. Can still be the same value
                conReach = combinations[sortedAngles[0]]
            else:
                conReach = combinations[sortedWidths[0]]
        else:
            # if all reach combinations have the same width differences. Check angle difference
            conReach = combinations[sortedAngles[0]]



    # if target ID is not in selected combination set connected reach to NaN
    if ID not in conReach:
        conReach = np.nan
    else:
        # select correct reach_id from combination
        CR = np.array(conReach)
        conReach = CR[CR != ID][0]

    return conReach

def get_connected_reach(row:gpd.GeoSeries, direction:str, df:gpd.GeoDataFrame, dfNode:gpd.GeoDataFrame, dfNodeConnection:pd.DataFrame)->Tuple[int|float, float, float]:
    """Extract down or upstream reach\n
    input:
    - row: selected row
    - direction: up or downstream direction (code does not work for upstream)
    - df: complete dataframe
    - dfNode: complete node dataframe
    - dfNodeConnection: junction Nodes with corresponding reach_id's\n
    output:
    - reach: connected reach. Can be reach_id or NaN
    - wn: wrong normal reach connection (does not work currently)
    - wo: wrong end reach
    """

    if (direction != 'up') & (direction != 'dn'):
        print('find_down_reach: wrong direction input')
    connCol   = f'rch_id_{direction}'
    numberCol = f'n_rch_{direction}'
    
    connectedReach    = row[connCol]


    junctionNode    = dfNodeConnection[dfNodeConnection.reach_id == row.reach_id]
    junctionNode    = junctionNode[junctionNode[direction] == True]
    if junctionNode.shape[0] != 0:
        junctionNode    = junctionNode.Node_id.values[0]
        junctionReaches = dfNodeConnection[dfNodeConnection.Node_id == junctionNode]
    else:
        junctionReaches = pd.DataFrame()

    wn,who = np.nan, np.nan

    if (junctionReaches.shape[0] > 2):
        reach = determine_connected_reach(row, df, dfNode, junctionReaches)
    elif row['n_rch_dn'] > 0:
            reach = connectedReach[0]
    else:
        reach = np.nan


    return reach, wn, who

def split_river_segments(df:gpd.GeoDataFrame,dfNode:gpd.GeoDataFrame, minReachLen:int = 4*12)->gpd.GeoDataFrame:
    """
    Split river segments according to width length relationship\n
    input:
    - df: complete dataframe
    - minReachLen: factor for minimal reach length\n
    output:
    - df: new data frame
    """
    df     = df.copy() # avoid chaining
    dfNode = dfNode.copy()
    
    newReaches  = [] # list for new combined reaches
    riverSegments = df.loc[df.include_flag == '0', 'river_segment'].unique()
    for rs in riverSegments:

        # select riversegment
        riverSegment = df[df.river_segment == rs].copy()

        # count number of occurances of reach_id's in reach_id and downstream connected columns
        dfConcat = pd.concat([riverSegment['reach_id'], riverSegment['dn_connected_reach']])
        dfConcat = dfConcat.value_counts(ascending = True).reset_index()
        endReach = dfConcat[dfConcat['count'] == 1]['index'] # reach with only one occurance is initial reach

        
        if endReach.shape[0] != 1:
            print(f'split_river_segments, Error in identifying end Reach. Should single value. River_segment: {rs}')
        
        # loop over all reaches in river segment
        shortReaches, short = [], False
        for i in range(riverSegment.shape[0]):
            
            if i == 0:
                rowID = endReach.values[0] # select starting reach_id
            else:
                rowID = newID # new ID
                
            # select current row and reach_id for new row
            row   = df[df.reach_id == rowID]
            df.loc[row.index, 'river_segment_order'] = i

            row   = row.iloc[0]
            newID = df[df.reach_id == rowID]['dn_connected_reach'].iloc[0]
            
            # skip ghost type reaches --> end reaches
            if row['type'] == 6:
                continue
            else:
                # check lenght requirement, or if previous reach was short
                if (row['reach_len'] < (minReachLen*row['node_mean_max_width'])) | (short == True):
                    short = True
                    
                    if len(shortReaches) > 0:
                        # select rows in new combined reach
                        dfShort      = df[df['reach_id'].isin(shortReaches)]
                        dfShortNodes = dfNode[dfNode['reach_id'].isin(shortReaches)]
                        checkLength  = dfShort['reach_len'].sum()

                        # check length for new reach based average max width
                        if (minReachLen * dfShortNodes['max_width'].mean()) > checkLength:
                            shortReaches.append(rowID) # if short add rowID to short reaches
                        else:
                            # if new reach meets requirement add to newReaches
                            newReaches.append(shortReaches) 

                            # check for current row id if short or not
                            if row['reach_len'] < (minReachLen*row['node_mean_max_width']):
                                shortReaches = [rowID]
                            else:
                                newReaches.append([rowID])
                                # reset values
                                short = False
                                shortReaches = []
                                

                    else: # add to short reaches if initial short reach
                        shortReaches.append(rowID)

                else: # if reach is not short keep as original
                    newReaches.append([rowID])



            if (i == (riverSegment.shape[0] -1)) & (short == True):
                # if final reach is reached but combined reach is short add to combined reaches
                newReaches.append(shortReaches)

    df.loc[df['river_segment_order'].isna(), 'river_segment_order'] = 0

    # create unstacked dataframe of combined reaches with combined_reach_id and reach_id
    dfNewReaches1 = pd.DataFrame(newReaches)
    dfNewReaches = dfNewReaches1.unstack().reset_index(name='reach_id')
    dfNewReaches = dfNewReaches.dropna().rename(columns={'level_0': 'reach_order',"level_1": "combined_reach_id"})
    dfNewReaches[['reach_id','combined_reach_id','reach_order' ]] = dfNewReaches[['reach_id','combined_reach_id','reach_order' ]].astype(int)

    # dfNR = dfT.copy()
    df     = df.merge(dfNewReaches, on = 'reach_id', how = 'left')
    df.loc[df['combined_reach_id'].isna(), 'reach_order'] = 0
    df.loc[df['combined_reach_id'].isna(), 'combined_reach_id'] = np.arange(df[df['combined_reach_id'].isna()].shape[0]) +\
                                                                            df['combined_reach_id'].max() + 1

    dfNode = dfNode.merge(dfNewReaches, on = 'reach_id', how = 'left')

    return df, dfNode

def find_dn_river_segment(df:'gpd.GeoDataFrame') -> 'gpd.GeoDataFrame':
    """
    create down stream river segment variable. For complete river topology.
    input:
    - df: total dataframe\n
    output:
    - dfN: updated dataframe
    """
    dfN = df.copy()

    dfSegment = df.loc[df.groupby('river_segment')['river_segment_order'].idxmax()]

    dfN = df.copy()
    for i, r in dfSegment.iterrows():
        segment   = r['river_segment']
        rchDN     = r['rch_id_dn_orig']
        reachInds = dfN[dfN['river_segment'] == segment].index
        
        dnSegment = np.nan
        if isinstance(rchDN, list):
            if len(rchDN) == 1:
                dnReach= rchDN[0]
            else:
                if np.isnan(r['dn_connected_reach']):
                    conRows = dfN[dfN['reach_id'].isin(rchDN)]
                    dnReach = conRows.loc[conRows['dist_out'].idxmin(), 'reach_id']
                else:
                    dnReach = r['dn_connected_reach']
            dnSegment = dfN[dfN['reach_id'] == dnReach].iloc[0]['river_segment']
            
        dfN.loc[reachInds, 'dn_river_segment'] = dnSegment
    
    return dfN

def split_files(df:'gpd.GeoDataFrame', threshold:'int') -> gpd.GeoDataFrame:
    """
    Split dataframe rows in different groups based on graph and interconnected reaches. If single group is to large create single group\n
    input:
    - df: dataframe
    - threshold: max size for file\n
    output:
    - df: updated dataframe
    """
    df = df.copy()
    edges = []
    for i, r in df.iterrows():
        DN = r['rch_id_dn_orig']
        if isinstance(DN, list):
            for dn in DN:
                edges.append([r['reach_id'], dn])
    
    G = nx.Graph()
    G.add_edges_from(edges)

    networks = list(nx.connected_components(G))
    dfN = pd.DataFrame(networks)
    dfN = dfN.unstack().reset_index(name='reach_id')[['level_1', 'reach_id']]
    dfN = dfN.sort_values('level_1').dropna().rename(columns={"level_1": "networkGraph"})

    dfNS = dfN.groupby('networkGraph', as_index = False).size()
    dfNS['cumSum'] = dfNS['size'].cumsum()

    dfNS['networkGroup'] = -1

    group, groupSize = [], 0
    for i, r in dfNS.iterrows():
        if r['size'] > threshold:
            dfNS.loc[i, 'networkGroup'] = dfNS['networkGroup'].max() + 1
        else:
            group.append(i)
            groupSize += r['size']
            if groupSize > threshold:
                dfNS.loc[group, 'networkGroup'] = dfNS['networkGroup'].max() + 1
                group, groupSize = [], 0

        if (i == dfNS.index[-1]) & (len(group) != 0):
            if groupSize > 1000:
                val = dfNS['networkGroup'].max() + 1
            else:
                val = dfNS['networkGroup'].max() 
            dfNS.loc[group, 'networkGroup'] = val

    dfN = dfN.merge(dfNS[['networkGraph', 'networkGroup']], how = 'left', on ='networkGraph')


    df = df.merge(dfN, how = 'left', on = 'reach_id')
    return df

def river_catchment_position(df):
    for N in df['networkGraph'].unique():    
        dfCatch = df[df['networkGraph'] == N]

        # create edges
        edges = []
        for i, r in dfCatch.iterrows():
            DN = r['rch_id_dn_orig']
            if isinstance(DN, list):
                for dn in DN:
                    edges.append([r['reach_id'], dn])
        
        G = nx.DiGraph()
        G.add_edges_from(edges)
        RG = G.reverse()
        
        startPoints = [node for node, out_degree in RG.out_degree() if out_degree == 0]

        dfSegments = dfCatch.groupby(['path_segs', 'pfaf2'], as_index = False).size()
        for i, rs in dfSegments.iterrows():
            dfSegment = dfCatch[(dfCatch['path_segs'] == rs['path_segs']) & 
                                (dfCatch['pfaf2'] == rs['pfaf2'])]
            segmentIndex = dfSegment.index
            segmentReach = dfSegment['reach_id'].iloc[0]
            
            upstreamReaches   = nx.ancestors(G, segmentReach)
            upstreamEndPoints = list(set(startPoints) & set(upstreamReaches))

            if len(upstreamEndPoints) > 0:
                rowID = dfCatch.loc[dfCatch['reach_id'].isin(upstreamEndPoints), 'dist_out'].idxmax()
                row   = dfCatch.loc[rowID]
                distVal   = row['dist_out']
                reachVal  = row['reach_id']

            else:
                print('Error no upstream end point detected')
                reachVal = np.nan
                distVal = dfSegment['dist_out'].iloc[0]
            df.loc[segmentIndex, 'max_dist_out'] = distVal
            df.loc[segmentIndex, 'up_reach_id']  = reachVal
    df['catchment_position'] = df['dist_out'] / df['max_dist_out']
    return df



def new_reach_definition(df, dfNode,min_reach_len_factor, directory, fileName, save = False):
    """Combine reaches based on width and angles between reaches"""

    maxBinSize = 5000

    df, dfNode       = filter_SWORD_input(df, dfNode)

    # select only included reaches for junction nodes
    dfNodeConnection = get_connection_nodes(df[df['include_flag'] == '0'])
    
    ############################################
    # Create downstream connection. Loop over all included reaches
    ############################################
    wrongNormalReaches, wrongHOReaches = [], []
    for i, r in df[df['include_flag'] == '0'].iterrows():
        IND = df[df.reach_id == r.reach_id].index[0]
        dn_reach, wrongNormal, wrongHO = get_connected_reach(r, 'dn', df, dfNode, dfNodeConnection)
        df.loc[IND, 'dn_connected_reach'] = dn_reach
        wrongNormalReaches.append(wrongNormal)
        wrongHOReaches.append(wrongHO)
    df['dn_connected_reach'] = df['dn_connected_reach'].astype(float)


    wrongNormalReaches =  dfNodeConnection[dfNodeConnection.Node_id.isin(wrongNormalReaches)].reach_id.values
    wrongHOReaches = np.array(wrongHOReaches)
    wrongHOReaches = wrongHOReaches[~np.isnan(wrongHOReaches)].astype(int)

    ############################################
    # add river segment ID to dataframe
    ############################################
    # create edges from down stream connections. 
        # Do not include reaches with no downstream connection, or with include_flag not 0
    edges = df[['reach_id', 'dn_connected_reach']]
    notIncludeReaches = df[df.include_flag != '0'].reach_id.values
    edges = edges[(~edges['dn_connected_reach'].isna()) & 
                    (~edges['reach_id'].isin(notIncludeReaches)) & 
                    (~edges['dn_connected_reach'].isin(notIncludeReaches))
                    ].values

    # create network from networkx to create all separate river segments
    G = nx.Graph()
    G.add_edges_from(edges)

    # create dataframe with river segments 
    riverSegments   = list(nx.connected_components(G)) # all individual connected segments
    dfriverSegments = pd.DataFrame(riverSegments)
    dfriverSegments = dfriverSegments.unstack().reset_index(name='reach_id')[['level_1', 'reach_id']]
    dfriverSegments = dfriverSegments.sort_values('level_1').dropna().rename(columns={"level_1": "river_segment"})
    # merge river segments with main dataframe
    df = df.merge(dfriverSegments, on = 'reach_id', how = 'left')
    # give each row a river_segment value. For reaches without value create new values
    df.loc[df['river_segment'].isna(), 'river_segment'] = np.arange(df[df['river_segment'].isna()].shape[0]) + df['river_segment'].max() + 1
    
    ############################################
    # add river segment ID to dataframe
    ############################################
    df, dfNode = split_river_segments(df, dfNode, min_reach_len_factor)

    # find downstream river segment
    df = find_dn_river_segment(df)

    ############################################
    # remove lakeflags
    ############################################
    def remove_lakeflag_ends(x):
        # Create a copy of the array for iteration
            arrIndices = np.arange(len(x))
            arrRemove = arrIndices.copy()
            for i in range(len(x)):
                if x[0] == 1:
                    x = np.delete(x, 0)
                    arrIndices = np.delete(arrIndices, 0)
                elif x[-1] == 1:
                    x = np.delete(x, -1)
                    arrIndices = np.delete(arrIndices, -1)
                else:
                    # Stop the loop if condition is not satisfied
                    break
            arrRemove = np.delete(arrRemove, arrIndices)
            return arrRemove

    # group by combined reaches and take mean of lakeflag id and check if combined reach contains lakeflag with value 1
    groupedLakeflag = df.groupby('combined_reach_id',as_index = False).agg(
                                    mean_lake_flag=('lakeflag', 'mean'),
                                    partial_lake =('lakeflag', lambda x: 1 in x.values))
    # select all combined reaches with only lakeflag 1 reaches
    lakes         = groupedLakeflag[(groupedLakeflag['mean_lake_flag'] % 1 == 0) & (groupedLakeflag['mean_lake_flag'] == 1)]['combined_reach_id'].values

    # combined reaches partially existing of lakeflags
    partial_lakes = groupedLakeflag[(groupedLakeflag['mean_lake_flag'] % 1 != 0) & (groupedLakeflag['partial_lake'] == True)]

    # select lakes at the endpoint of a combined reach
    removeReaches = []
    for cr in partial_lakes['combined_reach_id'].values:
        pl = df[df['combined_reach_id'] == cr].sort_values('reach_order')

        removeLake = remove_lakeflag_ends(pl['lakeflag'].values)
        removeReaches.extend(pl.iloc[removeLake].reach_id.values)

    # update include flag identified reaches
    removeReaches.extend(lakes)
    df = update_include_flag(df, removeReaches, '3')

    # remove lakes and partial lakes from combined_reach_id and assign new combined
    df.loc[df['include_flag'] == '3', 'combined_reach_id'] = np.arange(df[df['include_flag'] == '3'].shape[0]) + df['combined_reach_id'].max() + 1


    ############################################
    # remove short reaches?
    ############################################
    # remove reaches shorter than minimal requirement! 12*4*30
    df['combined_reach_len']   = df.groupby('combined_reach_id')['reach_len'].transform('sum')
    df['combined_reach_width'] = df.groupby('combined_reach_id')['width'].transform('mean')
    
    dfNode['combined_reach_max_width'] = dfNode.groupby('combined_reach_id')['max_width'].transform('mean')
    nmw = dfNode.groupby('combined_reach_id', as_index = False)['max_width'].mean()
    nmw = nmw.rename(columns = {'max_width':'combined_reach_max_width'})
    df = df.merge(nmw, on = 'combined_reach_id', how = 'left')
    
    shortReaches = df.loc[df['combined_reach_len'] < min_reach_len_factor*30, 'reach_id'].values
    df = update_include_flag(df, shortReaches, '4')

    ############################################
    # remove missing Up, Down, dn_connected reaches
    ############################################
    df = remove_updn_connection(df, True)


    ############################################
    # change values for reach_order after removal of reaches
    ############################################
    for cr in df.combined_reach_id.unique(): 
        rows = df[df.combined_reach_id == cr]
        rows = rows.sort_values('reach_order')
        df.loc[rows.index, 'reach_order'] = np.arange(rows.shape[0])

    ############################################
    # add include_flag to dfNode
    ############################################
    dfIncludeFlag = df[['reach_id', 'include_flag']].copy()
    dfNode = dfNode.merge(dfIncludeFlag, on = 'reach_id', how = 'left')
    

    ############################################
    # add catchment position
    ############################################
    df =  split_files(df, 5000)

    ############################################
    # add catchment position
    ############################################
    df = river_catchment_position(df)


    if save == True:
        dfSave = df.copy()
        dfSave[['rch_id_up', 'rch_id_dn', 
                'rch_id_up_orig', 'rch_id_dn_orig']] = dfSave[['rch_id_up', 'rch_id_dn',
                                                            'rch_id_up_orig', 'rch_id_dn_orig']].astype(str)
                
        for group in dfSave['networkGroup'].unique():
            dfGroupSave = dfSave[dfSave['networkGroup'] == group]
            dfNodeSave = dfNode[dfNode.reach_id.isin(dfGroupSave.reach_id.values)]
            groupString = f'{group}'
            if group < 10:
                groupString =  f'0{groupString}'
            
            dfGroupSave.to_file(directory + f'results/new_segments/vector/{fileName}_{groupString}_reach_new_segments.gpkg', driver = 'GPKG')
            dfNodeSave.to_file(directory + f'results/new_segments/node/{fileName}_{groupString}_node_new_segments.gpkg', driver = 'GPKG')

        
    return df, dfNode

