# Import whole packages
import geopandas as gpd
import pandas as pd
import shapely
import numpy as np
import networkx as nx



# Import functions from package
from shapely.geometry import Point
from scipy.optimize import minimize

from tqdm import tqdm

#Code
def spatial_self_join(DF):
    df = DF.copy()
    dfReduced = df[['reach_id', 'geometry']].copy()

    # add buffer for selfjoin with nearest neighbor
    dfReduced['geometry'] = dfReduced.buffer(300)

    # Self join of dataframe
    dfJoined = gpd.sjoin_nearest(dfReduced, dfReduced, 
                          max_distance = 10, distance_col = 'dist')
    
    # remove self match
    dfJoined = dfJoined[dfJoined.index != dfJoined.index_right]

    # group at left id --> all right id's are neighbors. And turn neighbors into list
    grouped        = dfJoined.groupby('reach_id_left')
    connected_list = {group : group_df.reach_id_right.tolist() for group, group_df in grouped}
    
    # create dataframe of neighbor list
    dfConnected = pd.DataFrame(list(connected_list.items()), columns=['Key', 'Value'])


    df = pd.merge(df, dfConnected, how = 'left', left_on = 'reach_id', right_on =  'Key')
    return df

def get_cycles(DF, DFNode):

    df = DF.copy()
    dfNode = DFNode.copy()
    
    # get edges based on grouped reach_ids
    grouped = dfNode.groupby('reach_id')
    # Get the indices of all rows that belong to each group
    edges_3 = [tuple(group_df.node_id.tolist()) for group, group_df in grouped if len(group_df.node_id.tolist()) > 2]
    edges_2 = [tuple(group_df.node_id.tolist()) for group, group_df in grouped if 
           (len(group_df.node_id.tolist()) > 1) &
           (len(group_df.node_id.tolist()) < 3)]
    
    for i in range(len(edges_3)):
        dfEdgeMany   = dfNode[dfNode.node_id.isin(edges_3[i])]
        mainReachId  = dfEdgeMany.reach_id.value_counts().keys()[0]
        mainReach    = df[df.reach_id == mainReachId].iloc[0]

        nodeRows = dfEdgeMany[dfEdgeMany.node_id.duplicated(keep = 'first') == False]
        dists = np.zeros(nodeRows.shape[0])
        for j in range(nodeRows.shape[0]):
            dists[j] = mainReach.geometry.project(nodeRows.iloc[j].geometry)

        sortedNodeIds = nodeRows.node_id.unique()[np.argsort(dists)]
        edgesMany     = [(sortedNodeIds[I], sortedNodeIds[I+1]) for I in range(len(sortedNodeIds) -1)]
    
        edges_2.extend(edgesMany)
    



    G = nx.Graph()
    G.add_edges_from(edges_2)
    cycles = nx.cycle_basis(G)
    
    cycleNodes = []
    cycleReaches = []
    for c in cycles:
        CNodes       = dfNode[dfNode.node_id.isin(c)]
        CNodesCounts = CNodes.reach_id.value_counts()

        
        cycleNodes.append(CNodes.node_id.values)
        cycleReaches.append(CNodesCounts[CNodesCounts > 1].keys())

    # find cycles with only 2 reaches
    res = list(set([ele for ele in edges_2 if edges_2.count(ele) > 1])) 

    for i, I in enumerate(res):
        rows = dfNode[dfNode.node_id.isin(I)]
        singleCycle = rows[rows.reach_id.duplicated()].reach_id.values
        cycleReaches.append(list(singleCycle))

    cycleDict = {}
    for i, cycle in enumerate(cycleReaches):
        for rid in cycle:
            rowIndex = df[df.reach_id == rid].index[0]

            if rowIndex in cycleDict:
                cycleDict[rowIndex].append(i)
            else:
                cycleDict[rowIndex] = [i]
    df['cycleID'] = df.index.map(cycleDict)
    df.loc[df['cycleID'].isna(), 'cycle'] = False
    df.loc[df['cycleID'].isna() == False , 'cycle'] = True    
        
    return df

def remove_wrong_connections(DF, wrongMid, DFNode):

    # add second midd group
    df = DF.copy()
    dfNode = DFNode.copy()
    for midComb in wrongMid:
      
        # remove mid's from vector dataframe
        for midIndex, mid in enumerate(midComb):
            r = df[df.reach_id == mid]
            row = r.iloc[0]

            for col in ['Mid', 'Mid1', 'Mid2']:
                if isinstance(row[col], list):
                    if (len(row[col]) == 1):
                        df.loc[r.index, col] = np.nan
                    else:

                        
                        Mid = row[col]

                        if midIndex == 0:
                            rem = midComb[1]
                        else:
                            rem = midComb[0]
                        

                        if rem in Mid:
                            Mid.remove(rem)

                            df.at[r.index[0], col] = list(Mid)


            # remove connection Node
            dropNodeId = dfNode[(dfNode.reach_id.isin(midComb)) & (dfNode.side == col) ].node_id
            if dropNodeId.shape[0] > 0:
                dfNode = dfNode[dfNode.node_id != dropNodeId.values[0]]
    

    return df, dfNode

def find_optimal_point(lines):
    # Objective function: sum of distances from a point (x, y) to all LineStrings
    def objective_function(coords):
        x, y = coords
        point = Point(x, y)

        # Sum the distances from the point to all lines
        total_distance = sum([point.distance(line) for line in lines])
        return total_distance
    

    # Initial guess for the point (e.g., at the centroid of the bounding box)
    initial_point = shapely.ops.nearest_points(lines[0], lines[1])[0].coords
    initial_guess = initial_point[0]

    # Use scipy's minimize to find the point that minimizes the sum of distances
    result = minimize(objective_function, initial_guess, method='Nelder-Mead')

    # The optimal point
    node = Point(result.x)
    return node

def nodes_for_reach(DF, reach:'int', S, M, M1, M2, E, dfNode, min_dist = 3e2):
    df = DF.copy()
    c = 0
    # print('start_nodes_for_reach: ', reach)
    sideName = ['Start', 'Mid','Mid1','Mid2', 'End']
    

    for sideIndex, side in enumerate([S, M,M1,M2, E]):
        

        reaches  = [reach]
        reachRow = df[df.reach_id == reach].iloc[0]
        
        if isinstance(side, list):

            
            reaches.extend(side)
            
            # find optimal point with global crs
            dfLines_glob = df[df.reach_id.isin(reaches)]
            lines_glob   = dfLines_glob.geometry.values
            node_glob    = find_optimal_point(lines_glob)
            
            # find optimal points with local crs (correct distance measuring)
            dfLines_loc  = dfLines_glob.copy()
            dfLines_loc  = dfLines_loc.to_crs(reachRow.localCRS)
            lines_loc    = dfLines_loc.geometry.values
            node_loc     = find_optimal_point(lines_loc)
            

            
            if dfNode.shape[0] > 0:
                NodeCheck = dfNode.copy()
                NodeCheck = NodeCheck.to_crs(reachRow.localCRS)
                nodeDistance = node_loc.distance(NodeCheck.geometry.values)
                nodeDistance = np.sort(nodeDistance)

                if nodeDistance[0] < min_dist: # add check if the distance is small if all are connected
                    continue;
            else:
                nodeDistance = [0]
        

            
            nodeId = dfNode.node_id.max() + 1
            if np.isnan(nodeId):
                nodeId = 0

            
            dfAdd = gpd.GeoDataFrame({'node_id' :[nodeId]*len(reaches),
                                      'reach_id': reaches,
                                      'geometry': [node_glob]*len(reaches),
                                      'distance': nodeDistance[0],
                                      'side'    :  sideName[sideIndex]  })
            dfNode = pd.concat([dfNode, dfAdd])           
            
            
    return dfNode

def check_mid_group(DF : gpd.GeoDataFrame, mid,reach, row, max_dist = 3e2):
    """
    input:
        DF: vector dataframe
        mid: list of Mid connected reaches
        row: dataseries row of main reach
        """
    midCheck = list(mid)
    midCheck.append(reach)

    dfMids = DF[DF.reach_id.isin(midCheck)].copy()
    dfMids = dfMids.to_crs(row.localCRS)
    P = []
    for I, R in dfMids.iterrows():
        P.append(shapely.ops.nearest_points(row.geometry, R.geometry)[1])
    
    
    if P[0].distance(P).max() > max_dist:
        diffMid = False
    else:
        diffMid = True
    return diffMid

def find_connected_side(DF, projection, fileName):
    """
    Identify all reaches that create a circular connection with connected reaches
    
    """
    df = DF.copy() # copy to avoid chaining dataframes
    
    # create self joined dataframe
    df = spatial_self_join(df)
    
    
    dfNode = gpd.GeoDataFrame(columns = ["node_id", 'reach_id', 'geometry', 'distance', 'side'])
    dfNode = dfNode.set_crs(projection)


    S, E, M, M1, M2 = [], [] ,[], [],[]
    MC, wrongMid = [], []
    # for i in tqdm(range(df.shape[0]), total = df.shape[0]):
    for i in range(df.shape[0]):
        S.append([])
        E.append([])
        M.append([])
        M1.append([])
        M2.append([])

        row = df.iloc[i]
        reaches = row['Value']
        if isinstance(reaches, list):
            for reach in reaches:
                # adjust code for multiple MIDS --> multiple midd connections

                dfLines = df[df.reach_id.isin([row.reach_id, reach])].copy()
                
                # use crs of main reach --> can be on the edge of local crs
                dfLines = dfLines.to_crs(row.localCRS) 

                rowCRS = dfLines[dfLines.reach_id == row.reach_id].iloc[0]
                targetCRS = dfLines[dfLines.reach_id == reach].iloc[0]
                lines = dfLines.geometry.values

                start = Point(rowCRS.geometry.coords[0])
                end   = Point(rowCRS.geometry.coords[-1])

                startDistance = start.distance(targetCRS.geometry)
                endDistance   = end.distance(targetCRS.geometry)
                
                # check connectionSide
                

                if startDistance < 4e2:
                    S[i].append(reach)

                if endDistance < 4e2:
                    E[i].append(reach)
                if (startDistance > 4e2) & (endDistance > 4e2):
                    # check for multiple mids and if they are different connections
                    if len(M[i]) > 0:
                        MCheck = check_mid_group(df, M[i],reach, row, 3e2)
                        
                        if MCheck == True:
                            M[i].append(reach)
                        elif len(M1[i]) > 0:
                            MCheck = check_mid_group(df, M1[i], reach, row,3e2)
                            if MCheck == True:
                                M1[i].append(reach)
                            elif len(M2[i]) > 0:
                                MCheck = check_mid_group(df, M2[i], reach, row,3e2)
                                if MCheck == True:
                                    M2[i].append(reach)
                                else:
                                    print(f'4 Mid connected reach!!!!!! ERROR {row.reach_id}')
                            else:
                                M2[i].append(reach)
                        else:
                            M1[i].append(reach)
                    else:
                        M[i].append(reach)

                    # check if Mids match --> wrong Mid. 2 centerlines cant be mid connected
                    midComb = list(np.sort([row.reach_id, reach]))
                    if midComb in MC:
                        wrongMid.append(midComb)
                    MC.append(midComb)   
                    
        
        ########################
        # Append Start
        ########################
        if len(S[i]) == 0:
            S[i] = np.nan
        

        ########################
        # Append Mid
        ########################
        #### MULTIGROUP!!!!!!!!!!!
        if len(M[i]) == 0:
            M[i] = np.nan
        if len(M1[i]) == 0:
            M1[i] = np.nan
        if len(M2[i]) == 0:
            M2[i] = np.nan
        

        ########################
        # Append End
        ########################
        if len(E[i]) == 0:
            E[i] = np.nan
        
        
            
    

        dfNode = nodes_for_reach(df, row.reach_id, S[i], M[i], M1[i], M2[i], E[i], dfNode, 3e2)   

    df['Start']  = S
    df['End']    = E
    df['Mid']    = M
    df['Mid1']   = M1
    df['Mid2']   = M2

    df, dfNode = remove_wrong_connections(df, wrongMid, dfNode)

    df = get_cycles(df, dfNode)
    df[['Key','Value','Start','End','Mid','Mid1','Mid2','cycleID']] = df[['Key','Value','Start','End','Mid','Mid1','Mid2','cycleID']].astype(str)
    
    # save cycle file with original geometries
    df.to_file(fileName +'_cycle.shp', driver='ESRI Shapefile', crs = projection)
    dfNode.to_file(fileName+'_connectionNodes.shp', driver='ESRI Shapefile', crs = projection)

    df = df[['reach_id', 'Key','Value','Start','End','Mid','Mid1','Mid2','cycleID','cycle']]

    return df, dfNode





