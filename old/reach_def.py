#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:20:14 2024

@author: niekcollotdescury
"""
import numpy as np
from shapely.geometry import Point, LineString, MultiLineString

import pandas as pd
import geopandas as gpd
import itertools


from line_functions import azimuth, localCRSReachLength

def find_starting_coords(lineList):

    # Identify starting and ending points for two linestrings
    ds = []
    ind = []
    L1 = list(lineList[0].coords)
    L2 = list(lineList[1].coords)
    for i in [0,-1]:
        Pi = Point(L1[i])
        for j in [0,-1]:
            Pj = Point(L2[j])
            ds.append(Pi.distance(Pj))
            ind.append([i,j])
    indices = ind[np.array(ds).argmin()]

    # Combine coordinates for new linestring
    if  indices == [-1,-1]:
        L2.reverse()
        L1.extend(L2)
    elif indices == [0 ,-1]:
        L1.reverse()
        L2.reverse()
        L1.extend(L2)
    elif indices == [-1, 0]:
        L1.extend(L2)
    elif indices == [0 , 0]:
        L1.reverse()
        L1.extend(L2)
        
    linest = LineString(L1)
    return linest

def multi_to_line(lineList):
    if len(lineList) == 2:

        outLine = find_starting_coords(lineList)
    elif len(lineList) == 3:

        intermediateLine = find_starting_coords(lineList)
        
        lineList = [intermediateLine, lineList[2]]
        outLine = find_starting_coords(lineList)

    return outLine


def too_short(DF, DF_node, projection, MW_ratio = 12, number_wave = 4):
    """Calculates if a reach length is too short based on the meander wavelength vs width of the river"""
    df = DF.copy()
    df_node = DF_node.copy()
    df_node = df_node[df_node.manual_add == 0]
    
    present = 'node_mwm' in df.columns
    if present == False:
        
        NMW = df_node.groupby('reach_id', as_index= False)['max_width'].mean()
        NMW = NMW.rename(columns={"max_width": "node_mwm"})
        df = pd.merge(df, NMW, on ='reach_id',how='left')

    
    df.loc[:, 'min_len'] = df.node_mwm * MW_ratio * number_wave
    df = localCRSReachLength(df, 'short', projection)
    df.loc[:, 'short'] = df.short < df.min_len
    return df




def find_connected_side(node1, node2):
    nodeS = node2.geometry.distance(node1.geometry.iloc[0]).min()
    nodeE = node2.geometry.distance(node1.geometry.iloc[-1]).min()
    if nodeS > nodeE:
        rowConnect = [-1, -2]
    else:
        rowConnect = [0,1]
    return rowConnect


# def azimuth(point1, point2):
#     '''azimuth between 2 shapely points (interval 0 - 360)'''
#     angle = np.arctan2(point2.x - point1.x, point2.y - point1.y)
#     return np.degrees(angle) if angle >= 0 else np.degrees(angle) + 360




def get_az(id, id_ud, df , df_node):
    row_ud = df[df.reach_id == id_ud]
    if row_ud.shape[0] > 0:
        row_ud = row_ud.iloc[0]

        
        nodes    = df_node[df_node.reach_id == id]
        nodes_ud = df_node[df_node.reach_id == row_ud.reach_id]
        
        if (nodes_ud.shape[0] < 2) | (nodes.shape[0] < 2):
            az = 1e10
        else: 
            nodesSide   = find_connected_side(nodes, nodes_ud)
            nodesSideUd = find_connected_side(nodes_ud, nodes)
            
            az_M  = azimuth(nodes.iloc[nodesSide[1]].geometry,nodes.iloc[nodesSide[0]].geometry)
            az_UD = azimuth(nodes_ud.iloc[nodesSideUd[0]].geometry,nodes_ud.iloc[nodesSideUd[1]].geometry)

            if (270 < az_M < 360) & (0 < az_UD < 90):
                az = (360 - az_M) + az_UD
            elif (270 < az_UD < 360) & (0 < az_M < 90):
                az = (360 - az_UD) + az_M
            else:
                az = np.abs(az_M - az_UD)


            
    else:
        az = 1e10
    return (az)


def check_manual_add(df_node,df_node_nm, id1, id2):
    #find the nodes for the selected reach and the up/down reach
    nodes    = df_node[df_node.reach_id == id1]
    nodesUd  = df_node[df_node.reach_id == id2]
    
    # find the nodes for the selected reach and the up/down 
    # reach for a dataframe where all manually added nodes are removed
    nodes_nm    = df_node_nm[df_node_nm.reach_id == id1]
    nodesUd_nm  = df_node_nm[df_node_nm.reach_id == id2]

    nodesSide   = find_connected_side(nodes, nodesUd)
    nodesSideUd = find_connected_side(nodesUd, nodes) 

    manAdd   = nodes.iloc[nodesSide[0]].node_id
    manAddUd = nodesUd.iloc[nodesSideUd[0]].node_id
    
    Add   = manAdd in nodes_nm.node_id.values
    AddUd = manAddUd in nodesUd_nm.node_id.values
    if (Add == False) | (AddUd == False):
        Add = True
    else:
        Add = False
    return Add



def confluence(id, df, ud_side):
    """Identify whether the down or upstream side is a confluence"""
    
    row = df[df.reach_id == id].iloc[0]
    if ud_side == 'up':
        rch_id = row.rch_id_up
    else:
        rch_id = row.rch_id_dn
    conf_type = 'No'
    conf_ids = []
    if len(rch_id) == 0:

        conf = False    
    elif len(rch_id) == 11:
        
        if ud_side == 'dn':
            row_dn = df[df.reach_id == int(rch_id)].iloc[0]
            if (row_dn.rch_id_up == None):
                conf = False

            else:
                if (len(row_dn.rch_id_up) <= 11):
                    conf = False

                else:

                    conf = True
                    conf_type = 'DN'
                    conf_ids = list(map(int, row_dn.rch_id_up.split(' ')))
                
        else:

            conf = False
    else:
        conf = True
        conf_type = 'Simple'
        conf_ids = list(map(int, rch_id.split(' ')))


    conf_ids.append(row.reach_id)
    return conf, conf_type, conf_ids



def connect_in_conf(id, conf_ids, df,df_node):

    comb = list(itertools.combinations(conf_ids,2))

    angles = np.zeros(len(comb))
    for i in range(len(comb)):
        angles[i] = get_az(comb[i][0], comb[i][1], df, df_node)

    connected = comb[angles.argmin()]
    connected_check = id in connected
    return connected_check


def identify_up_down(row, df,df_whole,df_node, updown, df_node_nm):
    """selected Row, vector dataframe, unfiltered vector dataframe, node dataframe, up or down indication """
    # max_angle = 60
    if updown =='up':
        # n_rch  = row.n_rch_up
        rch_id = row.rch_id_up
    else:
        # n_rch  = row.n_rch_dn
        rch_id = row.rch_id_dn   
    
    az_ret =np.nan
    if rch_id != None:

        if len(rch_id) > 11: # check if SWORD UD reach contains more than one UD reach
            # convert string of multiple reach id's into seperate ints
            ids_up = list(map(int, rch_id.split(' ')))
            rows_up = df[df.reach_id.isin(ids_up)] # Look for the id's in the dataframe

            # Update ids: remove ids that are not found in dataframe
            ids_up = rows_up.reach_id.values
            ids = []
            angles = []
            for ud in range(len(ids_up)):
                add      = check_manual_add(df_node, df_node_nm, row.reach_id, ids_up[ud])
                angle_ud = get_az(row.reach_id, ids_up[ud], df, df_node)
                if (add == False):
                    ids.append(ids_up[ud])
                    angles.append(angle_ud)
            
            
            if len(ids) == 0: # if no ud segments found return NaN
                rch_id = np.nan
            elif len(ids) == 1: # if one reach, select this reach if closest node is not mannually added
                conf, conf_type, conf_ids = confluence(row.reach_id, df_whole, updown)
                rch_id = ids[0]
                az_ret = angles[0]
                if (conf == True):
                    if conf_type == 'DN':
                        ##### NO DOWNSTREAM??
                        connected = connect_in_conf(row.reach_id, conf_ids, df,df_node)
                        if connected == False:
                            rch_id = np.nan
                            az_ret = np.nan                        
                    
                
            else:            
                angles   = np.array(angles)
                az_index = angles.argmin() 
                
                rch_id = ids[az_index]
                az_ret = angles[az_index]

        else:
            rch_id = df[df.reach_id == int(rch_id)]
            if rch_id.shape[0] == 0:
                rch_id = np.nan
            else:
                # CHECK PRESENTS OF UD SEGMENT
                # Check manual add
                # CHECK IF CONNECTION TO UD SEGMENT IS CONNECTED RIVER OR MAIN CHANNEL/TRIBUTARY
                rch_id = rch_id.iloc[0].reach_id
                
                add1 = check_manual_add(df_node, df_node_nm, row.reach_id, rch_id)
                angle_ud1 = get_az(row.reach_id, rch_id, df, df_node)

                conf, conf_type,conf_ids = confluence(row.reach_id, df_whole, updown)
                az_ret = angle_ud1
                if (add1 == True) :
                    rch_id = np.nan
                    az_ret = np.nan
                elif (conf == True):
                    if conf_type == 'DN':
                        connected = connect_in_conf(row.reach_id, conf_ids, df,df_node)
                        if connected == False:
                            rch_id = np.nan
                            az_ret = np.nan
                    elif conf_type == 'Simple':
                        print('\n\n\n\n\n\n\n\n ERROR: conf type should not be possible.\n\n\n\n\n\n\n\n\n')
                        
    else:
        rch_id = np.nan

    return rch_id, az_ret

def flatten_list(nested_list):
    flattened = []
    for item in nested_list:
        if isinstance(item, list):
            flattened.extend(flatten_list(item))  # Recursively flatten sub-lists
        else:
            flattened.append(item)  # If the item is not a list, just add it
    return flattened

def check_mid_connecting_reach(rid:'str', rch_id : 'str', df:'gpd.GeoDataFrame'):
    """
    input:
        rid: Target reach id
        rch_id: Up or down connected reach id
        df: Dataframe
    output:
        Connected reach id. If midd connected change to NaN
    """
    checkRow = df[df.reach_id == rch_id]

    midString = checkRow[['Mid', 'Mid1', 'Mid2']].values[0]
    for m in midString:
        if m != 'nan':
            mid = eval(m)
            if rid in mid:
                rch_id = np.nan


    return rch_id

def create_ud_id(DF,df_whole, DF_node, df_node_nm):

    df = DF.copy()
    df_node = DF_node.copy()

    df['up_segment'] = True
    df['dn_segment'] = True
    
    for i in range(df.shape[0]):

        r = df.iloc[i]     
        rch_id_up, azup = identify_up_down(r,df,df_whole, df_node, 'up', df_node_nm )
        rch_id_dn, azdn = identify_up_down(r,df,df_whole, df_node, 'dn', df_node_nm )

        if r.Value != 'nan':
            conVal = eval(r.Value)
            if rch_id_up not in conVal:
                rch_id_up = np.nan
                azup      = np.nan
            else:
                rch_id_up = check_mid_connecting_reach(r.reach_id, rch_id_up, df)
            if rch_id_dn not in conVal:
                rch_id_dn = np.nan
                azdn      = np.nan
            else:
                rch_id_dn = check_mid_connecting_reach(r.reach_id, rch_id_dn, df)
        else:
            rch_id_up, azup, rch_id_dn, azdn = np.nan, np.nan, np.nan, np.nan

        ind = df[df.reach_id == r.reach_id].index
        df.loc[ind, 'single_up'] = rch_id_up
        df.loc[ind, 'single_dn'] = rch_id_dn
        df.loc[ind, 'az_up'] = azup
        df.loc[ind, 'az_dn'] = azdn
    
        if np.isnan(rch_id_up):
            df.loc[ind, 'up_segment'] = False
        if np.isnan(rch_id_dn):
            df.loc[ind, 'dn_segment'] = False
    
    df['single_up'] = df.single_up.astype('Int64')
    df['single_dn'] = df.single_dn.astype('Int64')
    df['order'] = 0

    return df



def get_ud_id(row, df,dir, order):
    r = row
    for o in range(order):
        
        if dir == 'up':
            ud = r.single_up
        else:
            ud = r.single_dn
        if (ud is pd.NA) == False:
            ud_row = df[df.reach_id == ud].iloc[0]
            r = ud_row
            id = r.reach_id
        else:
            id = None
            break;

    return id

def connect_ud_geometries(DF, df_orig, order, projection):
    df = DF.copy()
    
    ids = df[(df.short) & ((df.up_segment) | (df.dn_segment))].reach_id.values
    
    for s in ids:
        rw = df[df.reach_id == s]
        r = rw.iloc[0]
            
        l = []

        for dir in ['up', 'dn']:
            
            id = get_ud_id(r,df,dir, order)
            

            
            if (id != None) & (id != r.reach_id):
                rupDF = df_orig[df_orig.reach_id == id]
                rupDF = localCRSReachLength(rupDF, 'RLength', projection)
                rup = rupDF.iloc[0]
                l.append(rup.geometry)
                
                
                df.loc[rw.index, 'LCheck'] += rup.RLength
                df.loc[rw.index, f'rch_ids_{dir}'] += f'{str(rup.reach_id)} '
            else:
                df.loc[rw.index, f'{dir}_segment'] = False

            # Always add row geometry after the check for upstream segment
            if dir == 'up':
                l.append(r.geometry)
    
        if len(l) > 1:
            new_geom = MultiLineString(l)
            
            if new_geom.geom_type == 'MultiLineString':
                   
                new_geom = multi_to_line(list(new_geom.geoms))
            df.loc[rw.index, 'geometry'] = new_geom
            df.loc[rw.index, 'order']    = order
            
    
    return df

def create_new_segments(DF, df_orig, DF_node, projection):
    df      = DF.copy()
    df_orig = DF.copy()
    df_node = DF_node.copy()

    # df.loc[:, 'LCheck'] = df.reach_len
    df = localCRSReachLength(df, 'LCheck', projection)
    df.loc[:, 'rch_ids_up'] = ''
    df.loc[:, 'rch_ids_dn'] = ''
    
    df = too_short(df,df_node, projection, 12, 4)
    
    short = df[df.short == True].shape[0]
    ud = df[(df.short == True) & (df.up_segment == False) & (df.dn_segment == False)].shape[0]
    order = 1

    while short != ud:
        df = connect_ud_geometries(df, df_orig, order, projection)
        df = too_short(df,df_node, projection, 12, 4)
    
        short = df[df.short == True].shape[0]
        ud    = df[(df.short == True) & (df.up_segment == False) & (df.dn_segment == False)].shape[0]
    
        if order == 35:
            short = 1
            ud    = 1
        
        order += 1
        
    df = localCRSReachLength(df, 'reach_len', projection)
    # df.loc[:, 'reach_len'] = df.geometry.length
    return df
