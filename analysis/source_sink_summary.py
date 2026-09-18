"""Summarize existing bend routing assignments; never runs a spatial overlay."""
import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from .source_sink import ROOT, CLASSES, join_assignments
import pyarrow.parquet as pq


def summarize(vectors_path, metrics_path, output_dir, class_map=None):
    OUT = Path(output_dir)
    OUT.mkdir(parents=True, exist_ok=True)
    cols = ['vector_id', 'source_file', 'file', 'reach_id', 'bendID', 'GMM_7_hard',
            'cp_height', 'cm_height', 'slope_inn', 'slope_out', 'slope_left', 'slope_right']
    extra_metrics = [c for c in ['bendHeight', 'lineSlope']
                     if c in pq.ParquetFile(vectors_path).schema_arrow.names]
    v = pd.read_parquet(vectors_path, columns=cols + extra_metrics)
    r = pd.read_parquet(metrics_path)
    d = join_assignments(v, r, class_map=class_map)
    audit = {'inputs': [str(vectors_path.resolve()), str(metrics_path.resolve())],
             'class_map': CLASSES if class_map is None else class_map,
             'joined_rows': len(d), 'unclassified_rows': int(d.confinement_class.isna().sum()),
             'unique_source_file_bendID': True, 'one_to_one_vector_id': True,
             'metric_definitions': {'cp_height': 'River center-point elevation (m), minimum within river bounds in sampled profiles.',
                                    'mean_lateral_slope': 'Arithmetic mean of existing slope_inn and slope_out (m/m); both required. These are raw lateral confinement slopes, not longitudinal slope or independent terrain relief.'}}
    audit['raw_metric_checks'] = {}
    for col in ['cp_height', 'slope_inn', 'slope_out'] + extra_metrics:
        a = pd.to_numeric(d[col], errors='coerce')
        audit['raw_metric_checks'][col] = {'missing': int(a.isna().sum()), 'negative': int(a.lt(0).sum()),
                                         'sentinels': int(a.isin([-9999, 99999]).sum()),
                                         'min': float(a.min()), 'max': float(a.max())}
        d[col] = a.replace([-9999, 99999, np.inf, -np.inf], np.nan)
    d['mean_lateral_slope'] = (d.slope_inn + d.slope_out) / 2
    metrics_to_summarize = ['cp_height', 'mean_lateral_slope', 'slope_inn', 'slope_out']
    metrics_to_test = ['cp_height', 'mean_lateral_slope']
    if 'bendHeight' in extra_metrics:
        metrics_to_summarize.append('bendHeight')
        metrics_to_test.append('bendHeight')
        audit['metric_definitions']['bendHeight'] = 'Mean DEM elevation along the bend (m).'
    if 'lineSlope' in extra_metrics:
        d['abs_lineSlope'] = d.lineSlope.abs()
        metrics_to_summarize.append('abs_lineSlope')
        metrics_to_test.append('abs_lineSlope')
        audit['metric_definitions']['abs_lineSlope'] = 'Magnitude of fitted centerline elevation gradient (m/m), shared by bends in the sampled combined reach. Not lateral confinement slope.'
    d = d[d.confinement_class.notna()].copy()
    d.to_parquet(OUT / 'bend_assignments_topography.parquet', index=False)
    summaries, tests, regions, reproduction = [], [], [], []
    for cohort in ['all_assigned', 'figure_zero_transitions']:
        base = d if cohort == 'all_assigned' else d[d.transitions.eq(0)]
        unc = base.confinement_class.eq('Unconfined')
        src = base.zone.eq('Source')
        groups = {
            'Unconfined_Source': base[unc & src],
            'OtherClasses_Source': base[~unc & src],
            'Unconfined_Global': base[unc],
            'Unconfined_NotSource': base[unc & ~src & base.zone.notna()],
            'Unconfined_BypassOrSink': base[unc & base.zone.isin(['Bypass', 'Sink'])],
        }
        for name, frame in groups.items():
            for metric in metrics_to_summarize:
                a = frame[metric].dropna()
                q = a.quantile([.25, .5, .75])
                summaries.append(dict(cohort=cohort, group=name, n=len(frame), metric=metric,
                                      n_valid=len(a), small_sample=len(a) < 300, n_missing=len(frame)-len(a),
                                      q25=q.loc[.25], median=q.loc[.5], q75=q.loc[.75],
                                      iqr=q.loc[.75]-q.loc[.25]))
        # Compare disjoint groups; the source subset and global total overlap.
        for ref in ['OtherClasses_Source', 'Unconfined_NotSource', 'Unconfined_BypassOrSink']:
            for metric in metrics_to_test:
                a = groups['Unconfined_Source'][metric].dropna()
                b = groups[ref][metric].dropna()
                t = mannwhitneyu(a, b, alternative='two-sided', method='asymptotic') if len(a) and len(b) else None
                tests.append(dict(cohort=cohort, metric=metric, reference=ref, n_source=len(a),
                                  n_reference=len(b), U=float(t.statistic) if t else np.nan, p=float(t.pvalue) if t else np.nan,
                                  p_display=('<1e-300' if t.pvalue < 1e-300 else f'{t.pvalue:.6g}') if t else 'not tested (empty group)',
                                  rank_biserial_source_minus_reference=2*t.statistic/(len(a)*len(b))-1 if t else np.nan))
        counts = groups['Unconfined_Source'].continent.value_counts()
        for cont, n in counts.items():
            regions.append(dict(cohort=cohort, continent=cont, n=int(n), pct=100*n/counts.sum()))
        selected = base[unc & base.zone.isin(['Source', 'Bypass', 'Sink'])]
        for zone, f in selected.groupby('zone'):
            reproduction.append(dict(cohort=cohort, zone=zone, n=len(f), total_n=len(selected),
                                     count_pct=100*len(f)/len(selected), length_m=f.length_m.sum(),
                                     length_pct=100*f.length_m.sum()/selected.length_m.sum()))
    pd.DataFrame(summaries).to_csv(OUT / 'group_summary.csv', index=False)
    pd.DataFrame(tests).to_csv(OUT / 'mann_whitney.csv', index=False)
    pd.DataFrame(regions).to_csv(OUT / 'regions.csv', index=False)
    pd.DataFrame(reproduction).to_csv(OUT / 'percentage_audit.csv', index=False)
    (OUT / 'audit.json').write_text(json.dumps(audit, indent=2))

    # Figure denominator and count/length shares for every confinement class.
    figure = d[d.figure_eligible].groupby(['confinement_class', 'zone']).agg(
        n=('vector_id', 'size'), sampled_length_m=('length_m', 'sum')).reset_index()
    for col in ['n', 'sampled_length_m']:
        figure[col + '_pct'] = 100 * figure[col] / figure.groupby('confinement_class')[col].transform('sum')
    figure.to_csv(OUT / 'figure_percentages.csv', index=False)
    return pd.DataFrame(summaries), pd.DataFrame(tests)


def plot_figure(output_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    frame = pd.read_csv(output_dir / 'figure_percentages.csv')
    table = frame.pivot(index='confinement_class', columns='zone', values='sampled_length_m_pct')
    table = table.reindex(index=['Confined', 'Asym partial', 'Sym partial', 'Unconfined'],
                          columns=['Sink', 'Bypass', 'Source']).fillna(0)
    ax = table.plot.bar(stacked=True, color=['#9dcfc6', '#fbf7bd', '#bebad7'], figsize=(8, 5))
    ax.set_ylabel('Sampled length (%)')
    ax.set_xlabel('Confinement class (zero-transition bends)')
    ax.legend(title='Zone', loc='upper left', bbox_to_anchor=(1.01, 1))
    for container in ax.containers:
        ax.bar_label(container, labels=[f'{v:.1f}%' if v else '' for v in container.datavalues], label_type='center')
    plt.xticks(rotation=25, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / 'source_sink.png', dpi=200)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vectors', type=Path, required=True)
    parser.add_argument('--metrics', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'results/source_sink')
    parser.add_argument('--plot', action='store_true', help='Recreate the sampled-length figure.')
    parser.add_argument('--class-map', type=Path, help='JSON mapping GMM labels to the four class names; defaults to historical manuscript labels.')
    args = parser.parse_args()
    class_map = {int(k): v for k, v in json.loads(args.class_map.read_text()).items()} if args.class_map else None
    summary, tests = summarize(args.vectors, args.metrics, args.output_dir, class_map)
    print(summary[summary.metric.isin(['cp_height', 'mean_lateral_slope'])].to_string(index=False))
    print(tests.to_string(index=False))
    if args.plot:
        plot_figure(args.output_dir)


if __name__ == '__main__':
    main()
