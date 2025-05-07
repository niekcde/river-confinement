#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:19:28 2024

@author: niekcollotdescury
"""
from shapely.geometry import Point

import matplotlib.pyplot as plt
import shapely
from line_functions import remove_man_add_section, remove_self_intersect
import pandas as pd

def remove_manual_add(DF, df_node):
    df = DF.copy()
    
    node_man = df_node[df_node.manual_add == 1].reach_id.unique()

    for id_val in node_man:
        # print(id_val)
        V_row = df[df.reach_id == id_val]
        # print(V_row.head())
        
        V = V_row.iloc[0]
        
        N = df_node[df_node.reach_id == id_val]
        NMA = N[N.manual_add == 1]
        NMN = N[N.manual_add == 0]
        
        N_S = N.node_id.min() in NMA.node_id.values
        N_E = N.node_id.max() in NMA.node_id.values
        NL = 'Not'
        if ((N_S == False) & (N_E  == False)) == False:
            if (N_S == True) & (N_E == True):
                if N.shape[0] == NMA.shape[0]:
                    df = df.drop(V_row.index, inplace = False)
                else:
                    node = N[N.node_id == NMN.node_id.min()].iloc[0].geometry
                    endP = Point(V.geometry.coords[0])
                    
                    NL = remove_man_add_section(V.geometry, node, endP)
    
                    node = N[N.node_id == NMN.node_id.max()].iloc[0].geometry
                    endP = Point(V.geometry.coords[-1])
                    
                    NL = remove_man_add_section(NL, node, endP)
            
            else:
                if (N_E == True):
    
                    node_id = NMN.node_id.max()
                    endP = Point(V.geometry.coords[-1])
                
                elif (N_S == True):
    
                    endP = Point(V.geometry.coords[0])
                    node_id = NMN.node_id.min()
    
            
                node = N[N.node_id == node_id].iloc[0].geometry
                NL = remove_man_add_section(V.geometry, node, endP)
            if NL != 'Not':
                df.loc[V_row.index, 'geometry'] = NL

    return df



def remove_man_add(DF, DF_node, plot = False):
    df_node = DF_node.copy()
    df = DF.copy()
    
    
    ManAddReaches = df_node[df_node.manual_add == 1].reach_id.unique()
    dfNMAR = df_node[df_node.reach_id.isin(ManAddReaches)]
    length = len(ManAddReaches)
    for i in range(length):
        ID = ManAddReaches[i]
        # print(ID)
        
        reach = dfNMAR.loc[dfNMAR.reach_id == ID, ['reach_id', 'node_id', 'manual_add', 'geometry']]
        reach = reach.sort_values('node_id')
        reach['value_grp'] = (reach.manual_add.diff(1) != 0).astype('int').cumsum()
        
        reachcon = pd.DataFrame({'ID' : reach.groupby('value_grp').reach_id.first(), 
                       'ManAdd': reach.groupby('value_grp').manual_add.first(),
                       'value_grp': reach.groupby('value_grp').value_grp.first(),
                       'value_grpF': reach.value_grp.min(),
                       'value_grpL': reach.value_grp.max(),
                       'totalSize': reach.shape[0],
                       'Consecutive' : reach.groupby('value_grp').size()}).reset_index(drop=True)
        reachcon['perc'] = reachcon.Consecutive / reachcon.totalSize
        
        reachCheckMiddle = reachcon[(reachcon.ManAdd == 1) & (reachcon.Consecutive > 3) 
                               & (reachcon.value_grp != reachcon.value_grpF ) & (reachcon.value_grp != reachcon.value_grpL)]
        
        reachCheckEnds = reachcon[(reachcon.ManAdd == 1) &
                                  ((reachcon.value_grp == reachcon.value_grpF ) | (reachcon.value_grp == reachcon.value_grpL))]
        
        # remove self intersecting parts of linestring
        Line = df[df.reach_id == ID].iloc[0].geometry
        Line = remove_self_intersect(Line)
        
        
        vector_row = df[df.reach_id == ID]
        
        # drop rows where line consists of only manual added nodes
        if reach[reach.manual_add == 0].shape[0] == 0:
            df = df.drop(vector_row.index)
            continue
    
        # remove starting and ending manually added points
        if reachCheckEnds.shape[0] > 0:
            FD = reachCheckEnds[reachCheckEnds.value_grp == reachCheckEnds.value_grpF]
            LD = reachCheckEnds[reachCheckEnds.value_grp == reachCheckEnds.value_grpL]
            
            # select starting and ending manually added nodes 
            startNodes = reach[(reach.value_grp == reachCheckEnds.value_grpF.values[0]) & (reach.manual_add == 1)]
            endNodes   = reach[(reach.value_grp == reachCheckEnds.value_grpL.values[0]) & (reach.manual_add == 1)]
            allNodes   = list(startNodes.node_id.values) + list(endNodes.node_id.values)
            
            
            if (startNodes.shape[0] == 1):
                startNode = reach.iloc[1].node_id
                allNodes += [reach.iloc[1].node_id]
            else:
                startNode = startNodes.node_id.max()
                
            if endNodes.shape[0] == 1:
                endNode = reach.iloc[-2].node_id
                allNodes += [reach.iloc[-2].node_id]
            else:
                endNode = endNodes.node_id.min()   
                
            
            if (FD.shape[0] > 0) & (LD.shape[0] > 0): # remove start and end

            
                newline1 = remove_man_add_section(Line, 
                                              reach[reach.node_id == startNode].iloc[0].geometry,
                                              reach[reach.node_id == reach.node_id.min()].iloc[0].geometry)
                newline = remove_man_add_section(newline1, 
                                              reach[reach.node_id == endNode].iloc[0].geometry,
                                              reach[reach.node_id == reach.node_id.max()].iloc[0].geometry)
            elif (FD.shape[0] > 0): # remove Start
                newline = remove_man_add_section(Line, 
                                              reach[reach.node_id == startNode].iloc[0].geometry,
                                              reach[reach.node_id == reach.node_id.min()].iloc[0].geometry)
            elif (LD.shape[0] > 0): # remove End
                
                newline = remove_man_add_section(Line, 
                                              reach[reach.node_id == endNode].iloc[0].geometry,
                                              reach[reach.node_id == reach.node_id.max()].iloc[0].geometry)
            
            df.loc[vector_row.index, 'geometry'] = newline
            drop_index = df_node[df_node.node_id.isin(allNodes)]
            df_node.drop(drop_index.index, inplace = True)
        
        # remove middle manually added points
        if reachCheckMiddle.shape[0] > 0:
            if reach.iloc[0].manual_add == 1:
                reach = reach[reach.value_grp != 1]

            if reach.iloc[-1].manual_add == 1:
                reach = reach[reach.value_grp != reach.value_grp.max()]

                
            nodes = reach[reach.value_grp.isin(reachCheckMiddle.value_grp.values)].sort_values('node_id').node_id.values
            node_min = nodes[0]
            node_max = nodes[-1]
            
            Lmin = remove_man_add_section(Line, 
                                          reach[reach.node_id == node_min].iloc[0].geometry,
                                          reach[reach.node_id == reach.node_id.max()].iloc[0].geometry)
            Lmax = remove_man_add_section(Line, 
                                          reach[reach.node_id == node_max].iloc[0].geometry,
                                          reach[reach.node_id == reach.node_id.min()].iloc[0].geometry)
            
            reachNodes = reach.node_id.sort_values().values
            if Lmin.length > Lmax.length:
                newline = Lmin
                
                newNodeIds = reachNodes[reachNodes <= node_min]
            else:
                newline = Lmax
                newNodeIds = reachNodes[reachNodes >= node_max]

                

            D = df_node.copy()
            remNodes = D[D.reach_id == ID]
            remNodes = remNodes[~remNodes.node_id.isin(newNodeIds)]
            df.loc[vector_row.index, 'geometry'] = newline
            
            df_node.drop(remNodes.node_id.index, inplace = True)
            # df_node = df_node[df_node.node_id.isin(newNodeIds)]
            
        if plot == True:
            TV = df[df.reach_id == ID]
            TN = df_node[df_node.reach_id == ID]
            
            f, ax = plt.subplots()
            plt.axis('equal')
            f.suptitle(f'{ID}')
            TV.plot(ax = ax, color = 'yellow', zorder = 10)
            TN.plot(ax = ax,marker='o', color='red', markersize = 10, zorder = 10)
            ax.plot(*Line.xy, zorder = 1, linewidth = 10 , color = 'black')
            plt.show()
            
    df_node = df_node[df_node.reach_id.isin(df.reach_id)]
    return df, df_node
