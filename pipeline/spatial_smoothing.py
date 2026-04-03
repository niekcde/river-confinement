#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse
import pickle
from datetime import datetime as dt
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm

from .paths import format_factor_token, load_project_paths
from .support import concat_nc_smooth_files


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

    grouped = multiIterator.groupby(['file', 'networkGroup'], as_index = False).agg({
        'networkGraph': list,
        'networkGroupSize': 'first',
    }).reset_index()
    return grouped.sort_values('networkGroupSize', ascending=False)[['file', 'networkGraph']].values


def bend_neighbor_graph(df, length_dict_file):
    nodes = []
    edges = []
    for network in tqdm(df['networkGraph'].unique(), position=0, leave=True):
        dfb = df[(df['networkGraph'] == network)]
        for _, row in dfb.iterrows():

            cri = int(row['combined_reach_id'])
            bri = row['bendRank']
            rows = dfb[(dfb['combined_reach_id'] == cri)]

            if rows.shape[0] == 1:
                upSide, dnSide = 'up', 'dn'
                upRankID, dnRankID = -1, 0
            elif bri == rows['bendRank'].min():
                upSide, dnSide = 'up', 'id'
                upRankID, dnRankID = -1, bri
            elif bri == rows['bendRank'].max():
                upSide, dnSide = 'id', 'dn'
                upRankID, dnRankID = bri - 2, 0
            else:
                upSide, dnSide = 'id', 'id'
                upRankID, dnRankID = bri, bri - 2

            upID = row[f'combined_reach_{upSide}']
            dnID = row[f'combined_reach_{dnSide}']

            neighbors = []
            lens = []
            if ~np.isnan(upID):
                if upID == cri:
                    upR = rows['bendRank'].iloc[upRankID]
                    lens.append(rows['bendLen'].iloc[upRankID])
                    neighbors.append([int(upID), upR])
                else:
                    upR = dfb.loc[(dfb['combined_reach_id'] == upID), ['bendLen', 'bendRank']]
                    if upR.shape[0] > 0:
                        lens.append(upR['bendLen'].iloc[upRankID])
                        upR = upR['bendRank'].iloc[upRankID]
                        neighbors.append([int(upID), upR])

            if ~np.isnan(dnID):
                if dnID == cri:
                    dnR = rows['bendRank'].iloc[dnRankID]
                    lens.append(rows['bendLen'].iloc[dnRankID])
                    neighbors.append([int(dnID), dnR])
                else:
                    dnR = dfb.loc[(dfb['combined_reach_id'] == dnID), ['bendLen', 'bendRank']]
                    if dnR.shape[0] > 0:
                        lens.append(dnR['bendLen'].iloc[dnRankID])
                        dnR = dnR['bendRank'].iloc[dnRankID]
                        neighbors.append([int(dnID), dnR])

            sid = f'{cri}_{bri}'
            nodes.append(sid)

            for i, neighbor in enumerate(neighbors):
                len1 = row['bendLen']
                len2 = lens[i]
                edges.append([sid, f'{neighbor[0]}_{neighbor[1]}', (len1 * 0.5) + (len2 * 0.5)])

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_weighted_edges_from(edges)
    length_dict = dict(nx.all_pairs_dijkstra_path_length(G, weight='weight'))

    with length_dict_file.open('wb') as handle:
        pickle.dump(length_dict, handle)

    return length_dict


def smooth_attributes(sid, attr, length_dict, df, max_dist=20000, max_neighbors=7):
    dist_map = length_dict[sid]
    dist_map[sid] = 0
    dist_map = dict(sorted(dist_map.items(), key=lambda item: item[1]))

    dist_keys = np.array(list(dist_map.keys()))
    dist_vals = np.array(list(dist_map.values()))

    mask = dist_vals < max_dist
    neighbors = dist_keys[mask][:max_neighbors]
    dists = dist_vals[mask][:max_neighbors]
    if len(neighbors) == 1:
        vals = df.loc[df['bendID'] == sid, attr].values[0]
        vals = np.concatenate([vals, np.array([0, 0, 0, 0])])
        return vals

    dfF = df.loc[df['bendID'].isin(neighbors), ['bendLen'] + ['bendID'] + attr]
    vals = dfF.set_index('bendID').loc[neighbors].reset_index().values

    dist_sigma = vals[:, 1][0] / 2
    weights = np.exp(- (dists ** 2) / (2 * dist_sigma ** 2))
    weights /= weights.sum()

    wmean = np.average(vals[:, 2:], weights=weights, axis=0)
    lenW = len(weights)
    wstd = np.sqrt(np.array(
        (weights @ (vals[:, 2:] - wmean) ** 2) / (((lenW - 1) / lenW) * np.sum(weights)),
        dtype=float,
    ))

    return np.concatenate([wmean, wstd])


def run_bend_smoothing(cont, df, *, single_smoothed_dir, cross_token, hf_token):
    print('Run bend Smoothing', cont)
    df = df.copy()
    if df.empty:
        return None

    df['bendRank'] = df.groupby('combined_reach_id')['bendDistOut'].rank(ascending=False).astype(int)
    df['bendID'] = df['combined_reach_id'].astype(int).astype(str) + '_' + df['bendRank'].astype(str)
    df = df.sort_values(['file', 'combined_reach_id', 'bendRank'])

    length_dict_file = single_smoothed_dir / f'length_dict_{cont}.pkl'
    ld = bend_neighbor_graph(df, length_dict_file)

    attr = [
        'slope_left_normalized', 'slope_right_normalized',
        'slope_out_normalized', 'slope_inn_normalized',
    ]
    attrS = [f'{a}_smooth' for a in attr] + [f'{a}_smoothSTD' for a in attr]
    for sid in tqdm(df['bendID'].values, position=0, leave=True):
        df.loc[df['bendID'] == sid, attrS] = smooth_attributes(sid, attr, ld, df, max_neighbors=5)

    print('Create smooth vars done')
    dfNC = df.to_xarray()
    ncFile = single_smoothed_dir / f'{cont}_{cross_token}_{hf_token}_smoothed.nc'
    if ncFile.exists():
        ncFile.unlink()

    dfNC.to_netcdf(ncFile)
    return ncFile


def _run_bend_smoothing_task(task):
    cont, df_cont, single_smoothed_dir, cross_token, hf_token = task
    return run_bend_smoothing(
        cont,
        df_cont,
        single_smoothed_dir=single_smoothed_dir,
        cross_token=cross_token,
        hf_token=hf_token,
    )


def _prepare_smoothing_dataframe(ds):
    dsSubset = ds[[
        'combined_reach_id', 'reach_id', 'file', 'bendLen', 'up_reach_id',
        'networkGraph', 'networkGroup', 'river_name',
        'dn_connected_reach', 'up_connected_reach', 'rch_id_dn', 'rch_id_dn_orig', 'rch_id_up', 'rch_id_up_orig',
        'combined_reach_up', 'combined_reach_dn',
        'ER_inn', 'ER_out', 'slope_inn', 'slope_out',
        'ER_left', 'ER_right', 'slope_left', 'slope_right',
        'apex', 'catchment_position', 'bendWidths', 'bendMaxWidths',
        'cp_height', 'cm_height', 'bendDistOut', 'max_dist_out', 'bendHeight',
        'bendSin', 'sin', 'ang',
        'conFactor', 'catchment_position', 'combined_reach_len'
    ]]
    df = dsSubset.to_dataframe().reset_index()

    df['ER_max_inout'] = df[['ER_inn', 'ER_out']].max(axis=1)
    df['slope_max_inout'] = df[['slope_inn', 'slope_out']].max(axis=1)

    df['slope_diff_bendSide'] = df['slope_inn'] - df['slope_out']
    df['slope_diff_bendSide_abs'] = abs(df['slope_inn'] - df['slope_out'])
    df['slope_diff_side'] = df['slope_left'] - df['slope_right']

    df['ER_diff_bendSide'] = df['ER_inn'] - df['ER_out']
    df['ER_diff_bendSide_abs'] = abs(df['ER_inn'] - df['ER_out'])
    df['ER_diff_side'] = df['ER_left'] - df['ER_right']

    df['nonDimAmp'] = df['apex'] / df['bendWidths']
    df['catch_height'] = df['bendHeight'] * df['catchment_position']

    df['slope_max'] = (df['cm_height'] - df['cp_height']) / (df['bendWidths'] * 0.5)
    df['slope_out_normalized'] = df['slope_out'] / df['slope_max']
    df['slope_inn_normalized'] = df['slope_inn'] / df['slope_max']
    df['slope_left_normalized'] = df['slope_left'] / df['slope_max']
    df['slope_right_normalized'] = df['slope_right'] / df['slope_max']

    df['bend_catch_pos'] = df['bendDistOut'] / df['max_dist_out']
    df['glob_reach_id'] = df.groupby(['file', 'combined_reach_id']).ngroup()
    return df


def run_spatial_smoothing(*, cross_factor=50, height_factor=2, config_path=None, workers=6, continents=None):
    paths = load_project_paths(config_path)
    paths.ensure_step7_dirs()

    cross_token = format_factor_token(cross_factor)
    hf_token = format_factor_token(height_factor)
    input_file = paths.single_values_dir / f'global_{cross_token}_{hf_token}_conf.nc'
    if input_file.exists() is False:
        raise FileNotFoundError(
            f"Step 7 input file not found at {input_file}. "
            "Run Step 6 aggregation first."
        )

    ds = xr.open_dataset(input_file)
    df = _prepare_smoothing_dataframe(ds)

    available_continents = [c for c in ['oc', 'as', 'sa', 'af', 'eu', 'na'] if c in set(df['file'].dropna())]
    if continents is None or len(continents) == 0:
        target_continents = available_continents
    else:
        target_continents = [c for c in continents if c in set(df['file'].dropna())]

    if len(target_continents) == 0:
        raise FileNotFoundError(
            f"No continent data were found in {input_file} for the requested Step 7 run."
        )

    tasks = [
        (cont, df[df['file'] == cont].copy(), paths.single_smoothed_dir, cross_token, hf_token)
        for cont in target_continents
    ]

    print('start Multi')
    dt1 = dt.now()
    if workers == 1 or len(tasks) == 1:
        continent_outputs = [_run_bend_smoothing_task(task) for task in tasks]
    else:
        with Pool(min(workers, len(tasks))) as pool:
            continent_outputs = pool.map(_run_bend_smoothing_task, tasks)

    print(f'Create nodes and edges finished: {dt.now() - dt1}')
    global_output = concat_nc_smooth_files(
        cross=cross_factor,
        hf=height_factor,
        config_path=config_path,
    )
    return {
        'continent_outputs': [output for output in continent_outputs if output is not None],
        'global_output': global_output,
    }


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Step 7 entrypoint: spatially smooth the aggregated confinement dataset.'
    )
    parser.add_argument(
        '--config',
        help='Path to config/paths.local.json. Defaults to config/paths.local.json when present.',
    )
    parser.add_argument(
        '--cross-factor',
        type=int,
        default=50,
        help='Cross-distance factor used in the current Step 7 run. Default: 50.',
    )
    parser.add_argument(
        '--height-factor',
        type=float,
        default=2,
        help='Height factor to smooth, for example 2, 1.5, or 0.5. Default: 2.',
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=6,
        help='Number of continent-level worker processes to use. Default: 6.',
    )
    parser.add_argument(
        '--continents',
        nargs='*',
        help='Optional subset of continent codes to smooth, for example oc af.',
    )
    return parser.parse_args(argv)


def main_cli(argv=None):
    args = parse_args(argv)
    run_spatial_smoothing(
        cross_factor=args.cross_factor,
        height_factor=args.height_factor,
        config_path=args.config,
        workers=args.workers,
        continents=args.continents,
    )


if __name__ == '__main__':
    main_cli()
