#%%
from datetime import datetime as dt
from tqdm import tqdm

import os
from glob import glob
import networkx as nx
import numpy as np
import pickle

# import packages
import xarray as xr
import geopandas as gpd
import numpy as np
import pandas as pd

from multiprocessing import Pool
from functools import partial

directory = '/scratch/6256481/'

import sys
sys.path.insert(0, directory + f'python/py_code/')
from support import concat_nc_smooth_files
#%%
def create_multiprocess_iterator(df):
    multiIterator = df.groupby(['file', 'networkGraph'], as_index = False)['index'].count()
    multiIterator = multiIterator.sort_values('index', ascending = False).reset_index()

    threshold = 5000
    multiIterator['networkGroup'] = -1

    for f in ['oc', 'af', 'sa', 'as', 'na', 'eu']:
        mi = multiIterator[multiIterator['file'] == f]
        group, groupSize = [], 0
        for i, r in mi.iterrows():
            if r['index'] > threshold:
                multiIterator.loc[i, 'networkGroup'] = multiIterator['networkGroup'].max() + 1
            else:
                group.append(i)
                groupSize += r['index']
                if groupSize > threshold:
                    multiIterator.loc[group, 'networkGroup'] = multiIterator['networkGroup'].max() + 1
                    group, groupSize = [], 0

            if (i == mi.index[-1]) & (len(group) != 0):
                if groupSize > threshold:
                    val = multiIterator['networkGroup'].max() + 1
                else:
                    val = multiIterator['networkGroup'].max() 
                multiIterator.loc[group, 'networkGroup'] = val
    multiIterator['networkGroupSize'] = multiIterator.groupby('networkGroup')['index'].transform('sum')
    multiIterator = multiIterator.sort_values('networkGroupSize', ascending = False)

    grouped = multiIterator.groupby(['file', 'networkGroup'], as_index = False).agg({'networkGraph': list,
                                                                                    'networkGroupSize':'first'}).reset_index()
    it = grouped.sort_values('networkGroupSize', ascending=False)[['file', 'networkGraph']].values
    return it

def bend_neighbor_graph(df):
    nodes = []
    edges = []
    for network in tqdm(df['networkGraph'].unique(), position=0, leave=True):
        dfb = df[(df['networkGraph'] == network)]
        for idx, row in dfb.iterrows():
            
            cri  = int(row['combined_reach_id'])
            bri  = row['bendRank']
            rows = dfb[(dfb['combined_reach_id'] == cri)]

            if rows.shape[0] == 1:
                upSide, dnSide = 'up', 'dn'
                upRankID, dnRankID = -1, 0
            elif bri == rows['bendRank'].min():
                upSide, dnSide = 'up', 'id'
                upRankID, dnRankID = -1, bri
            elif bri == rows['bendRank'].max():
                upSide, dnSide = 'id', 'dn'
                upRankID, dnRankID = bri -2, 0
            else:
                upSide, dnSide = 'id', 'id'
                upRankID, dnRankID = bri, bri -2

            upID = row[f'combined_reach_{upSide}']
            dnID = row[f'combined_reach_{dnSide}']

            neighbors = []
            lens      = []
            if ~np.isnan(upID):
                if upID == cri:
                    upR = rows['bendRank'].iloc[upRankID]
                    lens.append(rows['bendLen'].iloc[upRankID])
                    neighbors.append([int(upID),upR])
                else:
                    upR  = dfb.loc[(dfb['combined_reach_id'] == upID), ['bendLen', 'bendRank']]
                    if upR.shape[0] > 0:
                        lens.append(upR['bendLen'].iloc[upRankID])
                        
                        upR = upR['bendRank'].iloc[upRankID]
                        neighbors.append([int(upID),upR])

            if ~np.isnan(dnID):
                if dnID == cri:
                    dnR = rows['bendRank'].iloc[dnRankID]
                    lens.append(rows['bendLen'].iloc[dnRankID])
                    neighbors.append([int(dnID),dnR])
                else:
                    dnR  = dfb.loc[(dfb['combined_reach_id'] == dnID), ['bendLen', 'bendRank']]
                    if dnR.shape[0] > 0:
                        lens.append(dnR['bendLen'].iloc[dnRankID])

                        dnR = dnR['bendRank'].iloc[dnRankID]
                        neighbors.append([int(dnID),dnR])


            sid = f'{cri}_{bri}'
            nodes.append(sid)
            # if cri == 104485:
            #     print(sid, neighbors)

            for i, neighbor in enumerate(neighbors):
                len1 = row['bendLen']
                len2 = lens[i]
                edges.append([sid, f'{neighbor[0]}_{neighbor[1]}', (len1*0.5) + (len2*0.5)])

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_weighted_edges_from(edges)
    length_dict = dict(nx.all_pairs_dijkstra_path_length(G, weight='weight'))
    
    cont = df['file'].iloc[0]
    with open(directory + f'results/single_smoothed/length_dict_{cont}.pkl', 'wb') as f:
        pickle.dump(length_dict, f)

    return length_dict

# # 3. Smooth attributes using distance-weighted average (1 / path_length)
# 3. Smooth attributes using distance-weighted average (1 / path_length)
def smooth_attributes(sid, attr, length_dict, df, max_dist=20000, max_neighbors=7):
    # Get neighbors within a max path distance
    # dist_map = length_dict[sid]
    # dist_map[sid] = 0
    # # print(sid, dist_map)
    # # check that self is included
    # dist_items = [(k, v) for k, v in dist_map.items() if v <= max_dist]

    # dist_items = sorted(dist_items, key=lambda x: x[1])[:max_neighbors]
    # # print(dist_items)
    # if not dist_items:
    #     return df.loc[df['bendID'] == sid, attr].values[0]
    
    # lens = np.array([df.loc[df['bendID'] == k, 'bendLen'].values[0] for k, _ in dist_items])
    
    # dists1      = np.array([d for _, d in dist_items])  
    # dist_sigma1 = lens[0]/2
    # weights1 = np.exp(- (dists1 ** 2) / (2 * dist_sigma1 ** 2))
    # weights1 /= weights1.sum()

    
    # vals1 = []
    # for a in attr:
    #     vals1.append(np.array([df.loc[df['bendID'] == k, a].values[0] for k, _ in dist_items]))
    # dt10 = dt.now()
    # print(dt10-dt1)
    
    dist_map      = length_dict[sid]
    dist_map[sid] = 0
    dist_map      = dict(sorted(dist_map.items(), key=lambda item: item[1]))
    
    dist_keys = np.array(list(dist_map.keys()))
    dist_vals = np.array(list(dist_map.values()))

    mask      = dist_vals < max_dist
    neighbors = dist_keys[mask][:max_neighbors]
    dists     = dist_vals[mask][:max_neighbors]
    if len(neighbors) == 1:
        vals = df.loc[df['bendID'] == sid, attr].values[0]
        vals = np.concatenate([vals, np.array([0,0,0,0])])
        return vals

    dfF = df.loc[df['bendID'].isin(neighbors), ['bendLen'] +['bendID']+ attr]
    vals  = dfF.set_index('bendID').loc[neighbors].reset_index().values

    dist_sigma = vals[:, 1][0]/2
    weights = np.exp(- (dists ** 2) / (2 * dist_sigma ** 2))
    weights /= weights.sum()

    # weighted average
    wmean   = np.average(vals[:,2:],  weights=weights, axis = 0)
    lenW = len(weights)
    # Weighted std
    wstd = np.sqrt(np.array(( weights@(vals[:,2:]-wmean)**2) / (((lenW-1)/lenW)*np.sum(weights)) , dtype = float))

    return np.concatenate([wmean, wstd])


def run_bend_smoothing(cont, df):

    print('Run bend Smoothing', cont)
    df = df[df['file'] == cont].copy()
    # add rank to bends withing combinedReach
    df['bendRank'] = df.groupby('combined_reach_id')['bendDistOut'].rank(ascending=False).astype(int)
    df['bendID']   = df['combined_reach_id'].astype(int).astype(str) + '_' + df['bendRank'].astype(str)

    df = df.sort_values(['file', 'combined_reach_id', 'bendRank'])
    
    # ldFolder = directory + f'results/single_smoothed/length_dict_{cont}.pkl'
    # if os.path.exists(ldFolder):
    #     with open(ldFolder, 'rb') as f:
    #         ld = pickle.load(f)
    # else:
    ld = bend_neighbor_graph(df)


    attr = ['slope_left_normalized', 'slope_right_normalized',
                    'slope_out_normalized','slope_inn_normalized']
    attrS   = [f'{a}_smooth' for a in attr] + [f'{a}_smoothSTD' for a in attr]
    for sid in tqdm(df['bendID'].values, position=0, leave=True):
        df.loc[df['bendID'] == sid, attrS] = smooth_attributes(sid, attr, ld, df, max_neighbors=5)

    print('Create smooth vars done')
    dfNC = df.to_xarray()
    ncFile = directory + f'results/single_smoothed/{cont}_{cross}_{hf}_smoothed.nc'
    if os.path.exists(ncFile):
        os.remove(ncFile)

    dfNC.to_netcdf(ncFile)

#%%
cross = 50
hfList    = [2]
for hf in hfList:
    hf    = f'{hf}' if hf > 10 else f'0{hf}'

    f = directory + f'results/single_values/global_{cross}_{hf}_conf.nc'
    ds = xr.open_dataset(f)

    dsSubset = ds[[
        'combined_reach_id', 'reach_id', 'file', 'bendLen', 'up_reach_id', 
        'networkGraph', 'networkGroup', 'river_name',
        'dn_connected_reach', 'up_connected_reach', 'rch_id_dn', 'rch_id_dn_orig', 'rch_id_up', 'rch_id_up_orig',
        'combined_reach_up', 'combined_reach_dn',
        'ER_inn', 'ER_out', 'slope_inn', 'slope_out',
        'ER_left', 'ER_right', 'slope_left', 'slope_right',
        'apex', 'catchment_position', 'bendWidths', 'bendMaxWidths',
        'cp_height', 'cm_height','bendDistOut','max_dist_out','bendHeight',
        'bendSin','sin', 'ang', 
        'conFactor', 'catchment_position', 'combined_reach_len'
        ]]
    df = dsSubset.to_dataframe().reset_index()

    df['ER_max_inout']    = df[['ER_inn', 'ER_out']].max(axis = 1)
    df['slope_max_inout'] = df[['slope_inn', 'slope_out']].max(axis = 1)


    df['slope_diff_bendSide']     = df['slope_inn']     - df['slope_out']
    df['slope_diff_bendSide_abs'] = abs(df['slope_inn'] - df['slope_out'])
    df['slope_diff_side']         = df['slope_left']    - df['slope_right']


    df['ER_diff_bendSide']     = df['ER_inn']     - df['ER_out']
    df['ER_diff_bendSide_abs'] = abs(df['ER_inn'] - df['ER_out'])
    df['ER_diff_side']         = df['ER_left']    - df['ER_right']

    df['nonDimAmp'] = df['apex'] / df['bendWidths']

    df['catch_height'] = df['bendHeight'] * df['catchment_position']

    df['slope_max'] = (df['cm_height'] - df['cp_height']) / (df['bendWidths']*0.5)
    df['slope_out_normalized']   = df['slope_out']   / df['slope_max']
    df['slope_inn_normalized']   = df['slope_inn']   / df['slope_max']
    df['slope_left_normalized']  = df['slope_left']  / df['slope_max']
    df['slope_right_normalized'] = df['slope_right'] / df['slope_max']


    df['bend_catch_pos'] = df['bendDistOut'] / df['max_dist_out']
    df['glob_reach_id'] = df.groupby(['file', 'combined_reach_id']).ngroup()

    print('start Multi')
    dt1 = dt.now()
    if __name__ == '__main__':
        partial_func = partial(run_bend_smoothing, df=df)
        continents = ['oc', 'as', 'sa', 'af', 'eu', 'na']
        with Pool(6) as pool:
            pool.map(partial_func, continents)
            

    # multiResults = list(results)
    print(f'Create nodes and edges finished: {dt.now() - dt1}')

    concat_nc_smooth_files(directory, cross, hf)