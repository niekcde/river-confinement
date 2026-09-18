"""Prepare source/sink inputs from a new global NetCDF and saved GMM fit.

Separate metrics and geometry commands allow using existing analysis/GIS environments.
No model fitting or DEM resampling is performed.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

JOIN_KEYS = ['file', 'reach_id', 'combined_reach_id', 'apex', 'bendDistOut']


def fingerprint(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return {'path': str(path), 'bytes': path.stat().st_size,
            'mtime_ns': path.stat().st_mtime_ns, 'sha256': digest.hexdigest()}


def prepare_metrics(smoothed, assignments, output_dir):
    import xarray as xr
    fields = JOIN_KEYS + ['bendID', 'cp_height', 'cm_height', 'bendHeight',
                          'slope_inn', 'slope_out', 'slope_left', 'slope_right']
    with xr.open_dataset(smoothed) as ds:
        frame = ds[fields].to_dataframe().reset_index()
        attrs = {k: v.item() if hasattr(v, 'item') else v for k, v in ds.attrs.items()}
    labels = pd.read_parquet(assignments)
    keys = ['index', 'file', 'reach_id', 'combined_reach_id', 'bendID']
    merged = frame.merge(labels[keys + ['GMM_7_hard']], on=keys, how='outer',
                         validate='one_to_one', indicator=True)
    if merged['_merge'].eq('right_only').any() or merged.duplicated(JOIN_KEYS).any():
        raise ValueError('Assignments do not match global identifiers, or geometry keys are ambiguous.')
    audit = {'smoothed': fingerprint(smoothed), 'assignments': fingerprint(assignments),
             'global_rows': len(frame), 'classified_rows': len(labels),
             'all_assignment_identifiers_match': True, 'smoothing': attrs}
    output_dir.mkdir(parents=True, exist_ok=True)
    merged.drop(columns='_merge').rename(columns={'index': 'vector_id'}).to_parquet(
        output_dir / 'new_metrics.parquet', index=False)
    (output_dir / 'new_input_manifest.json').write_text(json.dumps(audit, indent=2))
    print(json.dumps(audit, indent=2), flush=True)


def prepare_geometry(metrics, geometry_dir, bends_dir, output_dir):
    import geopandas as gpd
    import pyarrow as pa
    import pyarrow.parquet as pq
    import shapely
    from .source_sink import split_dateline
    data = pd.read_parquet(metrics)
    files = sorted(geometry_dir.glob('??_??_50_02_conf.gpkg'))
    if len(files) != 43:
        raise ValueError(f'Expected 43 global HF2 geometry files; found {len(files)}')
    output = output_dir / 'vectors.parquet'
    writer = None
    seen = []
    records = []
    try:
        for path in files:
            code = path.name[:5]
            geo = gpd.read_file(path, columns=[k for k in JOIN_KEYS if k != 'file'])
            geo['file'] = code[:2]
            frame = geo.merge(data[data.file.eq(code[:2])], on=JOIN_KEYS, how='left',
                              validate='one_to_one', indicator=True)
            if not frame['_merge'].eq('both').all():
                raise ValueError(f'Geometry has unmatched bend identity: {path}')
            frame = frame.drop(columns='_merge')
            bend_path = bends_dir / f'{code}_50.parquet'
            bends = pd.read_parquet(bend_path, columns=[k for k in JOIN_KEYS if k != 'file'] + ['lineSlope'])
            bends['file'] = code[:2]
            frame = frame.merge(bends, on=JOIN_KEYS, how='left', validate='one_to_one', indicator=True)
            if not frame['_merge'].eq('both').all():
                raise ValueError(f'Bend table has unmatched identifiers: {bend_path}')
            frame = frame.drop(columns='_merge')
            seen.extend(frame.vector_id.tolist())
            frame['source_file'] = path.name
            frame['source_layer'] = path.stem
            geographic = gpd.GeoDataFrame(frame, geometry='geometry', crs=geo.crs).to_crs(4326)
            before_wkb = shapely.to_wkb(geographic.geometry.values)
            geographic.geometry = geographic.geometry.map(split_dateline)
            changed = before_wkb != shapely.to_wkb(geographic.geometry.values)
            projected = geographic.to_crs(8857)
            plain = pd.DataFrame(projected.drop(columns='geometry'))
            plain['geom_wkb'] = shapely.to_wkb(projected.geometry.values)
            table = pa.Table.from_pandas(plain, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(output, table.schema)
            writer.write_table(table)
            records.append({'geometry': fingerprint(path), 'bend_topography': fingerprint(bend_path),
                            'rows': len(plain), 'dateline_split_vector_ids':
                            geographic.loc[changed, 'vector_id'].astype(int).tolist()})
            print(code, len(plain), 'prepared', flush=True)
    finally:
        if writer is not None:
            writer.close()
    if len(seen) != len(set(seen)) or set(seen) != set(data.vector_id):
        raise ValueError('Prepared geometries do not cover global vector IDs exactly once.')
    (output_dir / 'geometry_manifest.json').write_text(json.dumps(
        {'rows': len(seen), 'crs': 'EPSG:8857', 'join_keys': JOIN_KEYS, 'inputs': records}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    metrics = commands.add_parser('metrics')
    metrics.add_argument('--smoothed', type=Path, required=True)
    metrics.add_argument('--assignments', type=Path, required=True)
    metrics.add_argument('--output-dir', type=Path, required=True)
    geometry = commands.add_parser('geometry')
    geometry.add_argument('--metrics', type=Path, required=True)
    geometry.add_argument('--geometry-dir', type=Path, required=True)
    geometry.add_argument('--bends-dir', type=Path, required=True)
    geometry.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'metrics':
        prepare_metrics(args.smoothed, args.assignments, args.output_dir)
    else:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        prepare_geometry(args.metrics, args.geometry_dir, args.bends_dir, args.output_dir)


if __name__ == '__main__':
    main()
