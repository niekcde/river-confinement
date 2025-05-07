#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 09:58:06 2024

@author: niekcollotdescury
"""
from support import get_ids
import numpy as np

def width_ratio_singleReach(DF_vector, DF_node):
    
    df_vector = DF_vector.copy()
    df_node = DF_node.copy()

    df_node = df_node[df_node.manual_add == 0]
    
    df_node['widthRatio'] = df_node.width / df_node.max_width
    df_node               = df_node[df_node.reach_id.isin(df_vector.reach_id.values)]
    df                    = df_node.groupby('reach_id', as_index = False)['widthRatio'].mean()

    df_v = df_vector.merge(df, how = 'left', on = 'reach_id')
    
    return df_v


def width_ratio(DF_vector, DF_node):
    """
    Wetted river width vs Total river width. Averaged per node
    input:
        DF_vector: Dataframe containing all reach information
        DF_node: Dataframe containing all node information
    
    output:
        - Vector dataframe with additional dimensionless widthRatio column
    """
    df_vector = DF_vector.copy()
    df_node = DF_node.copy()
    df_node = df_node[df_node.manual_add == 0]
    
    
    for i, r in df_vector.iterrows():
    # for i, r in tqdm(df_vector.iterrows(), total=df_vector.shape[0]):

        ids = get_ids(r) # use if you want to use entire reach
        # ids = [r.reach_id]
        if isinstance(r.rch_ids_up, str):
            ids.extend(list(map(int, r.rch_ids_up.split())))
        if isinstance(r.rch_ids_dn, str):
            ids.extend(list(map(int, r.rch_ids_dn.split())))

        Nodes      = df_node[df_node.reach_id.isin(ids)]
        widthRatio = Nodes.width / Nodes.max_width
        df_vector.loc[i, 'widthRatio'] = np.nanmean(widthRatio)

    
    return df_vector