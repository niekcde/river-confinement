"""Isolated Step 7 grid, legacy audit, and matched-bend comparison artifacts.

No production files are overwritten. Run --help for input and subset options.
"""
if __package__ in (None, ""):
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "pipeline"

import argparse
import itertools
import json
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score, davies_bouldin_score
from sklearn.mixture import GaussianMixture
from threadpoolctl import threadpool_limits

from .local_bend_smoothing import (ATTRIBUTES, OUTPUT_ATTRIBUTES, prepare_bends,
                                   build_topology, smooth_local, candidate_weights)
from .spatial_smoothing import _prepare_smoothing_dataframe, bend_neighbor_graph, smooth_attributes
from .clustering_confinement import prepare_confinement_clustering_dataframe

KEYS = ['file', 'networkGraph', 'bendID']


def load_input(path, continents=None):
    with xr.open_dataset(path) as ds:
        if continents:
            ds = ds.isel(index=np.flatnonzero(np.isin(ds['file'].values, continents)))
        df = _prepare_smoothing_dataframe(ds)
    if continents:
        df = df[df.file.isin(continents)].copy()
    if df.empty:
        raise ValueError('No input bends selected')
    return prepare_bends(df)


def legacy_frame(df, directory):
    """Call unchanged production functions, continent by continent."""
    frames = []
    audit = []
    for continent, frame in df.groupby('file', sort=False):
        frame = frame.copy()
        distances = bend_neighbor_graph(frame, directory / f'length_dict_{continent}.pkl')
        topology = build_topology(frame, 4)
        positions = {sid: pos for pos, sid in enumerate(frame.bendID)}
        for sid in frame.bendID:
            frame.loc[frame.bendID == sid, OUTPUT_ATTRIBUTES] = smooth_attributes(
                sid, ATTRIBUTES, distances, frame, max_neighbors=5)
        # Endpoint/interior examples spanning the selected networks.
        selected = pd.concat([frame.head(3), frame.tail(3), frame.iloc[[len(frame)//2]]]).drop_duplicates('bendID')
        for row in selected.itertuples():
            ids, actual, _, _, directions = candidate_weights(topology, positions[row.bendID], 4, .5, False)
            local = {frame.iloc[pos].bendID: (direction, distance)
                     for pos, direction, distance in zip(ids, directions, actual)}
            candidates = [(sid, dist) for sid, dist in sorted(distances[row.bendID].items(), key=lambda x:x[1])
                          if dist < 20000][:5]
            weights = np.exp(-0.5 * (np.array([d for _, d in candidates]) / (row.bendLen/2))**2)
            weights /= weights.sum()
            for (sid, distance), weight in zip(candidates, weights):
                direction, local_distance = local.get(sid, ('outside_directional_paths', np.nan))
                audit.append(dict(file=continent, networkGraph=row.networkGraph, focal=row.bendID,
                                  candidate=sid, distance=distance, weight=weight, is_self=sid==row.bendID,
                                  direction=direction, local_actual_distance=local_distance,
                                  path_distance_difference=local_distance-distance))
        frames.append(frame)
    pd.DataFrame(audit).to_csv(directory / 'legacy_neighbor_audit.csv', index=False)
    return pd.concat(frames)


def summaries(df, name):
    rows = []
    for continent, frame in [('global', df), *list(df.groupby('file'))]:
        for attr in OUTPUT_ATTRIBUTES:
            values = frame[attr]
            rows.append(dict(configuration=name, continent=continent, attribute=attr, n=len(values),
                             n_finite=int(np.isfinite(values).sum()), mean=values.mean(), std=values.std(ddof=0),
                             q05=values.quantile(.05), median=values.median(), q95=values.quantile(.95)))
    return rows


def slope_comparison(left, right, a, b):
    joined = left.set_index(KEYS)[OUTPUT_ATTRIBUTES + ['bendLen']].join(
        right.set_index(KEYS)[OUTPUT_ATTRIBUTES], how='inner', lsuffix='_a', rsuffix='_b', validate='one_to_one')
    # Paired comparisons use a common bend population; quartiles come from input lengths.
    groups = [('global', joined)] + [(str(c), g) for c,g in joined.groupby(level='file')]
    bins = pd.qcut(joined.bendLen, 4, duplicates='drop')
    groups += [(f'length:{interval}', joined.loc[index]) for interval,index in bins.groupby(bins, observed=True).groups.items()]
    rows = []
    for group, frame in groups:
        for attr in OUTPUT_ATTRIBUTES:
            x, y = frame[attr+'_a'], frame[attr+'_b']
            valid = np.isfinite(x) & np.isfinite(y)
            delta = y[valid] - x[valid]
            rows.append(dict(reference=a, configuration=b, group=group, attribute=attr,
                             matched=len(frame), valid=int(valid.sum()), mean_change=delta.mean(),
                             mean_absolute_change=delta.abs().mean(), median_absolute_change=delta.abs().median(),
                             p95_absolute_change=delta.abs().quantile(.95),
                             rmse=np.sqrt(np.mean(delta**2)), max_absolute_change=delta.abs().max()))
    return rows


def clustering_grid(frames, output, clusters, seeds, sample_size, models):
    """Refit the existing eight-feature/PCA preparation with fixed model settings.

    Eligibility is reported before restricting to the intersection. All fits and
    scoring samples then use identical bend IDs and ordering across settings.
    """
    feature_cols = ATTRIBUTES + [a+'_smooth' for a in ATTRIBUTES]
    eligible = {}
    counts = []
    for name, frame in frames.items():
        indexed = frame.set_index(KEYS).sort_index()
        mask = np.isfinite(indexed[feature_cols].to_numpy(dtype=float)).all(axis=1)
        eligible[name] = indexed.loc[mask]
        counts.append(dict(configuration=name, total=len(frame), eligible=int(mask.sum())))
    common = None
    for frame in eligible.values():
        common = frame.index if common is None else common.intersection(frame.index)
    common = common.sort_values()
    for row in counts:
        row['common_fit_population'] = len(common)
    pd.DataFrame(counts).to_csv(output/'clustering_eligibility.csv', index=False)
    if len(common) < 3:
        raise ValueError('Too few common valid bends for clustering')
    invalid = [k for k in clusters if not 1 < k < len(common)]
    if invalid:
        raise ValueError(f'Invalid cluster counts for {len(common)} common bends: {invalid}')
    scores, comparisons, proportions = [], [], []
    assignments = {}
    for name, frame in eligible.items():
        prepared, cols = prepare_confinement_clustering_dataframe(frame.loc[common].reset_index())
        x = prepared[cols].to_numpy()
        labels_file = prepared[KEYS].copy()
        for model, k, seed in itertools.product(models, clusters, seeds):
            if model == 'kmeans':
                estimator = KMeans(n_clusters=k, n_init=10, random_state=seed)
            else:
                estimator = GaussianMixture(n_components=k, covariance_type='tied', n_init=5,
                                            max_iter=200, init_params='kmeans', random_state=seed)
            with threadpool_limits(limits=1):
                labels = estimator.fit_predict(x)
            tag = f'{model}_k{k}_seed{seed}'
            labels_file[tag] = labels
            assignments[(name, tag)] = labels
            sample = np.random.default_rng(seed).choice(len(x), min(sample_size,len(x)), replace=False)
            sampled_labels = labels[sample]
            valid_score = 1 < len(np.unique(sampled_labels)) < len(sample)
            scores.append(dict(configuration=name, model=model, k=k, seed=seed,
                               silhouette=silhouette_score(x[sample],sampled_labels) if valid_score else np.nan,
                               davies_bouldin=davies_bouldin_score(x[sample],sampled_labels) if valid_score else np.nan,
                               converged=getattr(estimator,'converged_',True)))
        labels_file.to_csv(output/f'{name}_cluster_assignments.csv', index=False)
    for a,b in itertools.combinations(frames,2):
        for model,k,seed in itertools.product(models,clusters,seeds):
            tag=f'{model}_k{k}_seed{seed}'
            la,lb=assignments[(a,tag)],assignments[(b,tag)]
            table=np.zeros((k,k),dtype=int)
            np.add.at(table,(la,lb),1)
            rows,cols=linear_sum_assignment(-table)
            mapping=np.empty(k,dtype=int)
            mapping[cols]=rows
            aligned=mapping[lb]
            for group, mask in [('global',np.ones(len(common),dtype=bool)),
                                *[(c,np.asarray(common.get_level_values('file')==c)) for c in common.get_level_values('file').unique()]]:
                comparisons.append(dict(reference=a,configuration=b,model=model,k=k,seed=seed,group=group,
                                        ari=adjusted_rand_score(la[mask],lb[mask]),
                                        changed_fraction=float(np.mean(la[mask]!=aligned[mask]))))
                for label in range(k):
                    pa,pb=np.mean(la[mask]==label),np.mean(aligned[mask]==label)
                    proportions.append(dict(reference=a,configuration=b,model=model,k=k,seed=seed,group=group,
                                            cluster=label,reference_proportion=pa,proportion=pb,change_pp=100*(pb-pa)))
    for filename, records in [('clustering_scores',scores),('clustering_comparison',comparisons),('cluster_proportions',proportions)]:
        pd.DataFrame(records).to_csv(output/f'{filename}.csv',index=False)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--continents',nargs='+')
    p.add_argument('--neighbors',type=int,nargs='+',default=[2,3,4])
    p.add_argument('--alphas',type=float,nargs='+',default=[.5,.75,1.])
    p.add_argument('--baseline',type=Path,help='Existing baseline NetCDF; compared on matched bend IDs')
    p.add_argument('--run-legacy',action='store_true',help='Compute unchanged legacy baseline on the selected input')
    p.add_argument('--max-network-bends',type=int,help='Diagnostic subset only: exclude larger networks')
    p.add_argument('--network-limit',type=int,help='Diagnostic subset only: retain this many complete networks')
    p.add_argument('--no-floor-control',nargs=2,type=float,metavar=('N','ALPHA'))
    p.add_argument('--clusters',nargs='+',type=int,help='Enable fixed-setting clustering comparisons for these k values')
    p.add_argument('--seeds',nargs='+',type=int,default=[20,43,50])
    p.add_argument('--models',nargs='+',choices=['kmeans','gmm'],default=['kmeans','gmm'])
    p.add_argument('--sample-size',type=int,default=4000)
    args=p.parse_args(argv)
    if args.baseline and args.run_legacy:
        p.error('Choose an existing baseline or --run-legacy')
    if any(n<1 for n in args.neighbors) or any(not np.isfinite(a) or a<=0 for a in args.alphas):
        p.error('Neighbors and alphas must be positive')
    args.output_dir.mkdir(parents=True,exist_ok=False)
    (args.output_dir/'manifest.json').write_text(json.dumps({**vars(args),'python':platform.python_version(),
        'status':'started','clustering_settings':{'kmeans_n_init':10,'gmm_n_init':5,'gmm_covariance':'tied','gmm_max_iter':200}},default=str,indent=2))
    start=time.perf_counter()
    df=load_input(args.input,args.continents)
    if args.max_network_bends or args.network_limit:
        sizes=df.groupby(['file','networkGraph']).size().sort_values(ascending=False)
        if args.max_network_bends:
            sizes=sizes[sizes<=args.max_network_bends]
        if args.network_limit:
            sizes=sizes.head(args.network_limit)
        mask=pd.MultiIndex.from_frame(df[['file','networkGraph']]).isin(sizes.index)
        df=df.loc[mask].copy()
    if df.empty:
        raise ValueError('No bends remain after network filtering')
    timings=[dict(stage='load_and_select',configuration='shared',seconds=time.perf_counter()-start,bends=len(df))]
    df.groupby(['file','networkGraph']).size().rename('bends').to_csv(args.output_dir/'network_sizes.csv')
    frames={}
    if args.run_legacy:
        start=time.perf_counter()
        frames['legacy']=legacy_frame(df,args.output_dir)
        timings.append(dict(stage='legacy_graph_and_average',configuration='legacy',seconds=time.perf_counter()-start,bends=len(df)))
    elif args.baseline:
        with xr.open_dataset(args.baseline) as ds:
            if args.continents:
                ds = ds.isel(index=np.flatnonzero(np.isin(ds['file'].values,args.continents)))
            baseline=ds.to_dataframe().reset_index()
        if 'bendID' not in baseline:
            baseline=prepare_bends(baseline)
        baseline=baseline.set_index(KEYS)
        frames['legacy']=baseline.loc[baseline.index.intersection(pd.MultiIndex.from_frame(df[KEYS]))].reset_index()
        if frames['legacy'].empty:
            raise ValueError('No baseline bend IDs match the input')
    start=time.perf_counter()
    depth=max(args.neighbors)
    if args.no_floor_control:
        n,a=args.no_floor_control
        if not n.is_integer() or n<1 or not np.isfinite(a) or a<=0:
            p.error('Control requires a positive integer N and positive alpha')
        depth=max(depth,int(n))
    topology=build_topology(df,depth)
    timings.append(dict(stage='adjacency_and_cached_walks',configuration='shared',seconds=time.perf_counter()-start,bends=len(df)))
    (args.output_dir/'topology_diagnostics.json').write_text(json.dumps(topology.diagnostics,indent=2))
    settings=[(n,a,True) for n,a in itertools.product(args.neighbors,args.alphas)]
    if args.no_floor_control:
        settings.append((int(args.no_floor_control[0]),args.no_floor_control[1],False))
    for n,a,floor in settings:
        name=f'n{n}_a{a:g}_floor{int(floor)}'
        start=time.perf_counter()
        frames[name]=smooth_local(df,topology,neighbors=n,alpha=a,length_floor=floor)
        timings.append(dict(stage='average',configuration=name,seconds=time.perf_counter()-start,bends=len(df)))
        ds=frames[name].to_xarray()
        ds.attrs.update(smoothing_method='directional_gaussian',neighbors_per_direction=n,
                        smoothing_alpha=a,smoothing_length_floor=int(floor),source_input=str(args.input))
        start=time.perf_counter()
        ds.to_netcdf(args.output_dir/f'{name}.nc')
        timings.append(dict(stage='write',configuration=name,seconds=time.perf_counter()-start,bends=len(df)))
        # Retain only comparison columns across the grid, not repeated metadata.
        frames[name] = frames[name][KEYS + ['bendLen'] + ATTRIBUTES + OUTPUT_ATTRIBUTES].copy()
        audit = []
        for focal in sorted({0, len(df)//2, len(df)-1}):
            ids, actual, effective, weights, directions = candidate_weights(topology,focal,n,a,floor)
            for pos, d, de, weight, direction in zip(ids,actual,effective,weights,directions):
                audit.append(dict(focal=df.iloc[focal].bendID, file=df.iloc[focal]['file'],
                                  networkGraph=df.iloc[focal].networkGraph, candidate=df.iloc[pos].bendID,
                                  direction=direction, actual_distance=d, effective_distance=de, weight=weight))
        pd.DataFrame(audit).to_csv(args.output_dir/f'{name}_neighbor_audit.csv',index=False)
        print(name, 'complete',flush=True)
    if args.run_legacy:
        frames['legacy'].to_xarray().to_netcdf(args.output_dir/'legacy.nc')
    pd.DataFrame(timings).to_csv(args.output_dir/'timings.csv',index=False)
    pd.DataFrame([row for name,frame in frames.items() for row in summaries(frame,name)]).to_csv(args.output_dir/'slope_summaries.csv',index=False)
    pd.DataFrame([row for a,b in itertools.combinations(frames,2) for row in slope_comparison(frames[a],frames[b],a,b)]).to_csv(args.output_dir/'slope_comparisons.csv',index=False)
    if args.clusters:
        clustering_grid(frames,args.output_dir,args.clusters,args.seeds,args.sample_size,args.models)
    manifest=json.loads((args.output_dir/'manifest.json').read_text())
    manifest['status']='complete'
    manifest['bends']=len(df)
    (args.output_dir/'manifest.json').write_text(json.dumps(manifest,indent=2))


if __name__=='__main__':
    main()
