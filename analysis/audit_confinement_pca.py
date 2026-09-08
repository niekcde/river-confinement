"""Read-only audit of archived Fh=2 input using the production PCA function."""
import json
import sqlite3
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import sklearn
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import clustering_confinement as production

SOURCE = Path('/Volumes/PhD/confinement/results/single_values/global_50_02_smoothed.nc')
EXPORT = SOURCE.parent.parent / 'dfClusterSlope.csv'
OUT = ROOT / 'results/pca_confinement_audit'
FEATURES = [f'slope_{s}_normalized' for s in ('out', 'right', 'left', 'inn')]
FEATURES += [f'{c}_smooth' for c in FEATURES.copy()]
IDS = ['file', 'reach_id', 'bendID']

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with xr.open_dataset(SOURCE) as ds:
        metadata = {'source': str(SOURCE), 'source_bytes': SOURCE.stat().st_size,
                    'source_attrs': ds.attrs, 'sklearn_version': sklearn.__version__}
        df = ds[FEATURES + IDS + ['cm_height', 'cp_height', 'bendWidths', 'conFactor'] +
                [f'slope_{s}' for s in ('out', 'right', 'left', 'inn')]].to_dataframe().reset_index()
    metadata['input_rows'] = len(df)
    metadata['input_continents'] = df.file.value_counts().to_dict()
    inferred_hf = (df.cm_height - df.cp_height) / (df.bendWidths * df.conFactor)
    metadata['height_factor_from_values'] = {'count': int(inferred_hf.notna().sum()),
        'min': float(inferred_hf.min()), 'max': float(inferred_hf.max())}
    metadata['notebook_normalization_max_abs_difference'] = {}
    maximum = (df.cm_height - df.cp_height) / (df.bendWidths * .5)
    for side in ('out', 'right', 'left', 'inn'):
        diff = df[f'slope_{side}_normalized'] - df[f'slope_{side}'] / maximum
        metadata['notebook_normalization_max_abs_difference'][side] = float(diff.abs().max())
    # Capture the actual fitted estimator without changing its implementation,
    # arguments, input features, filtering, scaling, or returned scores.
    fitted = []
    original = production.PCA
    def capture(*args, **kwargs):
        pca = original(*args, **kwargs)
        fitted.append(pca)
        return pca
    with patch.object(production, 'PCA', capture):
        clustered, pc_cols = production.prepare_confinement_clustering_dataframe(df)
    pca = fitted[0]
    metadata['retained_rows'] = len(clustered)
    metadata['retained_continents'] = clustered.file.value_counts().to_dict()
    metadata['solver'] = pca._fit_svd_solver
    metadata['explained_variance_ratio'] = pca.explained_variance_ratio_.tolist()
    metadata['cumulative_first_three'] = float(pca.explained_variance_ratio_[:3].sum())
    loadings = pd.DataFrame(pca.components_[:3].T, index=FEATURES, columns=['PC1','PC2','PC3'])
    loadings.to_csv(OUT / 'loadings_pc1_pc3.csv', index_label='variable')
    print(loadings.to_string(), flush=True)
    print(json.dumps(metadata, indent=2), flush=True)
    # Match by bend identity, since export rows may differ in order or coverage.
    export = pd.read_csv(EXPORT, usecols=IDS + FEATURES + pc_cols)
    metadata['export_rows'] = len(export)
    metadata['export_continents'] = export.file.value_counts().to_dict()
    metadata['duplicate_ids_input'] = int(clustered.duplicated(IDS).sum())
    metadata['duplicate_ids_export'] = int(export.duplicated(IDS).sum())
    joined = clustered[IDS + FEATURES + pc_cols].merge(export, on=IDS, how='outer',
        suffixes=('_input','_export'), indicator=True, validate='one_to_many')
    metadata['identity_match_counts'] = joined['_merge'].value_counts().to_dict()
    both = joined[joined['_merge'] == 'both']
    metadata['matched_max_abs_difference'] = {c: float((both[c+'_input']-both[c+'_export']).abs().max()) for c in FEATURES + pc_cols}
    metadata['matched_mean_abs_difference'] = {c: float((both[c+'_input']-both[c+'_export']).abs().mean()) for c in FEATURES + pc_cols}
    # Refit exported feature matrix through the same production function as well.
    with patch.object(production, 'PCA', capture):
        production.prepare_confinement_clustering_dataframe(export.drop_duplicates(IDS))
    export_pca = fitted[-1]
    metadata['export_explained_variance_ratio'] = export_pca.explained_variance_ratio_.tolist()
    pd.DataFrame(export_pca.components_[:3].T, index=FEATURES, columns=['PC1','PC2','PC3']).to_csv(OUT / 'export_loadings_pc1_pc3.csv', index_label='variable')
    (OUT / 'audit.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata, indent=2), flush=True)

def check_submission():
    source = SOURCE.parent.parent / 'submission/v1/global_clusters_bend.gpkg'
    with sqlite3.connect(f'file:{source}?mode=ro&immutable=1', uri=True) as connection:
        submitted = pd.read_sql_query('select ' + ','.join(IDS + FEATURES) +
            ' from global_clusters_bend where file is not null', connection)
    with xr.open_dataset(SOURCE) as ds:
        original = ds[IDS + FEATURES].to_dataframe().dropna(subset=FEATURES)
    matched = submitted.merge(original, on=IDS, how='outer',
        suffixes=('_submission', '_input'), indicator=True, validate='one_to_one')
    report = {'submission_file': str(source), 'submission_complete_rows': len(submitted),
        'matches': matched['_merge'].value_counts().to_dict(),
        'max_abs_feature_difference': {f: float((matched[f+'_submission']-
            matched[f+'_input']).abs().max()) for f in FEATURES}}
    (OUT / 'submission_comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)

if __name__ == '__main__':
    main()
    check_submission()
