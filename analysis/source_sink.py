"""Martin & Lamb routing-zone overlay, preserving the manuscript sampling method.

Run with ``python -m analysis.source_sink --help``. No work runs on import.
The saved assignment tables should be reused for manuscript comparisons.
"""
from __future__ import annotations

import argparse
import json
import math
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import shapely
from shapely.geometry import box
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
ZONES = {0: 'Ocean', 1: 'Sink', 2: 'Bypass', 3: 'Source', 4: 'No data'}
CLASSES = {0: 'Unconfined', 2: 'Confined', 5: 'Sym partial',
           1: 'Asym partial', 3: 'Asym partial', 4: 'Asym partial', 6: 'Asym partial'}
METHOD = 'manuscript_legacy_endpoint_v1'


def split_dateline(geometry):
    """Split geographic bend lines at ±180 before projecting to Equal Earth.

    Preserve original vertices; interpolate only the dateline intersection.
    Returns the original geometry unchanged when no longitude jump occurs.
    """
    from shapely.geometry import LineString, MultiLineString
    from shapely.ops import split
    from shapely.affinity import translate
    original_parts = line_parts(geometry)
    if not any(np.any(np.abs(np.diff(np.asarray(line.coords)[:, 0])) > 180)
               for line in original_parts):
        return geometry
    result = []
    for line in original_parts:
        coords = np.asarray(line.coords).copy()
        coords[:, 0] = np.rad2deg(np.unwrap(np.deg2rad(coords[:, 0])))
        pieces = [LineString(coords)]
        lo, hi = coords[:, 0].min(), coords[:, 0].max()
        for k in range(math.floor((lo - 180) / 360), math.ceil((hi - 180) / 360) + 1):
            boundary = 180 + 360 * k
            if lo < boundary < hi:
                cutter = LineString([(boundary, -90), (boundary, 90)])
                pieces = [part for piece in pieces for part in line_parts(split(piece, cutter))]
        for piece in pieces:
            mid = (piece.bounds[0] + piece.bounds[2]) / 2
            shift = -360 * math.floor((mid + 180) / 360)
            result.append(translate(piece, xoff=shift))
    return MultiLineString(result)


def require_unique_ids(frame):
    if 'vector_id' not in frame or frame.vector_id.isna().any() or not frame.vector_id.is_unique:
        raise ValueError('vector_id must be present, non-null and unique; never regenerate IDs for saved metrics.')


def load_vectors(path, geometry=False):
    """Read historical plain WKB or native GeoParquet without changing row IDs."""
    names = pq.ParquetFile(path).schema_arrow.names
    columns = names if geometry else [c for c in names if c not in ('geom_wkb', 'geometry')]
    frame = pd.read_parquet(path, columns=columns)
    require_unique_ids(frame)
    if geometry:
        geom_col = 'geom_wkb' if 'geom_wkb' in frame else 'geometry'
        if geom_col not in frame:
            raise ValueError('Input vectors need WKB geometry.')
        # Historical vectors.parquet was explicitly written in EPSG:8857.
        crs = 'EPSG:8857'
        meta = pq.ParquetFile(path).schema_arrow.metadata or {}
        if b'geo' in meta:
            geo = json.loads(meta[b'geo'])
            crs = geo['columns'][geom_col].get('crs', 'OGC:CRS84')
            if crs is None:
                raise ValueError('GeoParquet has an unknown CRS.')
        geoms = shapely.from_wkb(frame.pop(geom_col).map(bytes).to_numpy())
        frame = gpd.GeoDataFrame(frame, geometry=geoms, crs=crs).to_crs('EPSG:8857')
        if frame.geometry.isna().any() or frame.geometry.is_empty.any():
            raise ValueError('Empty/missing vector geometry; repair inputs without renumbering IDs.')
    return frame


def join_assignments(vectors, metrics, class_map=None):
    class_map = CLASSES if class_map is None else class_map
    unknown = set(vectors.GMM_7_hard.dropna()) - set(class_map)
    if unknown:
        raise ValueError(f'Unmapped GMM labels: {sorted(unknown)}')
    require_unique_ids(vectors)
    require_unique_ids(metrics)
    if not set(metrics.dominant_cat.dropna()).issubset(ZONES):
        raise ValueError('Unexpected routing category; expected 0, 1, 2, 3, 4.')
    joined = vectors.merge(metrics, on='vector_id', how='outer', validate='one_to_one', indicator=True)
    if not joined['_merge'].eq('both').all():
        raise ValueError('Vector and metric IDs do not match completely. Use the original paired tables.')
    if joined.duplicated(['source_file', 'bendID']).any():
        raise ValueError('Duplicate bend identities (source_file, bendID).')
    joined = joined.drop(columns='_merge')
    joined['confinement_class'] = joined.GMM_7_hard.map(class_map)
    joined['zone'] = joined.dominant_cat.map(ZONES)
    joined['continent'] = joined.source_file.str[:2]
    joined['figure_eligible'] = (joined.transitions.eq(0) & joined.confinement_class.notna()
                                 & joined.zone.isin(['Sink', 'Bypass', 'Source']))
    return joined


def raster_paths(path):
    """Read mask tiles directly from the supplied ZIP, or an extracted directory."""
    path = Path(path).resolve()
    if path.suffix.lower() == '.zip':
        with zipfile.ZipFile(path) as archive:
            names = sorted(n for n in archive.namelist()
                           if Path(n).match('mask_strat_*.tif'))
        paths = [f'/vsizip/{path}/{name}' for name in names]
    else:
        paths = [str(p) for p in sorted(path.rglob('mask_strat_*.tif'))]
    if not paths:
        raise ValueError(f'No mask_strat_*.tif tiles in {path}')
    return paths


def sample_categories(line, src, step):
    """Historical endpoint sampling: each sample carries length/(n-1) weight."""
    if line.length == 0:
        return np.array([], dtype=float), 0.
    n = max(2, int(math.ceil(line.length / step)) + 1)
    points = shapely.line_interpolate_point(line, np.linspace(0, line.length, n))
    coords = shapely.get_coordinates(points)
    samples = src.sample(coords)
    values = np.asarray(samples if isinstance(samples, np.ndarray) else list(samples)).reshape(-1)
    return values, line.length / (n - 1)


class CachedRasterSampler:
    """Same nearest-cell values as rasterio.sample, with one tile read into RAM."""
    def __init__(self, dataset):
        self.values = dataset.read(1)
        self.transform = dataset.transform
        self.nodata = dataset.nodata
        self.height, self.width = self.values.shape

    def sample(self, coords):
        from rasterio.transform import rowcol
        rows, cols = rowcol(self.transform, coords[:, 0], coords[:, 1])
        rows, cols = np.asarray(rows, dtype=np.int64), np.asarray(cols, dtype=np.int64)
        valid = (rows >= 0) & (rows < self.height) & (cols >= 0) & (cols < self.width)
        result = np.full(len(coords), self.nodata if self.nodata is not None else 0,
                         dtype=self.values.dtype)
        result[valid] = self.values[rows[valid], cols[valid]]
        return result[:, None]


def sample_metrics(values, segment_length, nodata):
    valid = np.isfinite(values)
    if nodata is not None:
        valid &= values != nodata
    cats = {}
    for value in values[valid]:
        if value not in ZONES:
            raise ValueError(f'Unexpected raster category {value}')
        key = int(value)
        cats[key] = cats.get(key, 0.) + segment_length
    # Legacy transitions skip nodata, including a gap between valid samples.
    sequence = values[valid]
    return cats, int(np.count_nonzero(np.diff(sequence))), int(valid.sum()), len(values)


def line_parts(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type == 'LineString':
        return [geometry]
    if hasattr(geometry, 'geoms'):
        return [line for part in geometry.geoms for line in line_parts(part)]
    return []


def overlay(vectors, tiles, step=50.):
    """Reproduce legacy metrics, retaining even vectors with no valid samples.

    Tile order is deterministic. Tied dominant categories use first encountered
    category, as in the original script. Transitions are counted per clipped part.
    """
    import rasterio

    if not np.isfinite(step) or step <= 0:
        raise ValueError('Sample step must be finite and positive.')
    require_unique_ids(vectors)
    if vectors.crs is None or vectors.crs.to_epsg() != 8857:
        raise ValueError('Overlay vectors must use EPSG:8857.')
    accum = {int(v): {'cats': {}, 'transitions': 0, 'length': 0., 'coverage': 0.}
             for v in vectors.vector_id}
    mapping = []
    index = vectors.sindex
    for tile_id, tile in enumerate(tqdm(tiles, desc='Routing tiles')):
        with rasterio.open(tile) as src:
            if src.crs is None or src.crs.to_epsg() != 8857 or src.count != 1:
                raise ValueError(f'Expected single-band EPSG:8857 raster: {tile}')
            bounds = box(*src.bounds)
            sampler = CachedRasterSampler(src)
            for pos in sorted(index.query(bounds, predicate='intersects')):
                row = vectors.iloc[pos]
                vid = int(row.vector_id)
                mapping.append({'vector_id': vid, 'tile_id': tile_id, 'tile': tile})
                state = accum[vid]
                tile_length = tile_coverage = 0.
                for line in line_parts(row.geometry.intersection(bounds)):
                    vals, seg = sample_categories(line, sampler, step)
                    cats, trans, nvalid, ntotal = sample_metrics(vals, seg, src.nodata)
                    length = nvalid * seg
                    tile_length += length
                    # Original coverage expression simplifies to valid sampled length.
                    tile_coverage += (nvalid / ntotal) * ntotal * seg if ntotal else 0.
                    state['transitions'] += trans
                    for cat, weight in cats.items():
                        state['cats'][cat] = state['cats'].get(cat, 0.) + weight
                state['length'] += tile_length
                coverage = min(1., tile_coverage / tile_length) if tile_length else 0.
                state['coverage'] += coverage * tile_length
    records = []
    for vid, state in accum.items():
        cats, length = state['cats'], state['length']
        proportions = np.asarray(list(cats.values())) / length if length else np.array([])
        record = {'vector_id': vid, 'length_m': length,
                  'coverage_pct': min(1., state['coverage'] / length) if length else 0.,
                  'transitions': state['transitions'],
                  'transitions_per_km': state['transitions'] / (length / 1000.) if length else 0.,
                  'entropy': float(-np.sum(proportions * np.log(proportions))),
                  'dominant_cat': max(cats, key=cats.get) if cats else np.nan}
        record.update({f'cat_{cat}_pct': cats.get(cat, 0.) / length if length else 0. for cat in ZONES})
        records.append(record)
    return pd.DataFrame(records), pd.DataFrame(mapping, columns=['vector_id', 'tile_id', 'tile'])


def prepare_vectors(folder, pattern, output):
    """Prepare a NEW paired run. Existing manuscript vector IDs must be preserved."""
    parts = []
    for path in sorted(folder.glob(pattern)):
        for layer in gpd.list_layers(path)['name']:
            frame = gpd.read_file(path, layer=layer)
            if frame.crs is None:
                raise ValueError(f'Missing CRS: {path}:{layer}')
            frame = frame[frame.geometry.geom_type.isin(['LineString', 'MultiLineString'])].copy()
            if frame.empty:
                continue
            frame = frame.to_crs('EPSG:4326')
            frame.geometry = frame.geometry.map(split_dateline)
            frame = frame.to_crs('EPSG:8857')
            frame['source_file'], frame['source_layer'] = path.name, layer
            parts.append(frame)
    if not parts:
        raise ValueError(f'No line features matching {folder / pattern}')
    frame = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs='EPSG:8857')
    if frame.duplicated(['source_file', 'bendID']).any():
        raise ValueError('Duplicate bend identities in input layers.')
    frame['vector_id'] = np.arange(len(frame), dtype=np.int64)
    frame.to_parquet(output, index=False)
    return len(frame)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    prep = commands.add_parser('prepare', help='Prepare vectors for a new overlay; IDs differ from historical ordering.')
    prep.add_argument('--vector-dir', type=Path, required=True)
    prep.add_argument('--glob', default='??_??_bend_cluster.gpkg')
    prep.add_argument('--output', type=Path, required=True)
    run = commands.add_parser('overlay', help='Explicitly run raster sampling (potentially expensive).')
    run.add_argument('--vectors', type=Path, required=True)
    run.add_argument('--rasters', type=Path, default=ROOT / 'data/source_sink/mask_strat_241022.zip')
    run.add_argument('--output-dir', type=Path, default=ROOT / 'results/source_sink_overlay')
    run.add_argument('--sample-step', type=float, default=50.)
    args = parser.parse_args()
    if args.command == 'prepare':
        if args.output.exists():
            parser.error('Output already exists; use a new path to protect paired vector IDs.')
        args.output.parent.mkdir(parents=True, exist_ok=True)
        print(f'Prepared {prepare_vectors(args.vector_dir, args.glob, args.output):,} bends')
    else:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        dest = args.output_dir / 'vector_raster_metrics.parquet'
        if any((args.output_dir / f).exists() for f in
               ['vector_raster_metrics.parquet', 'vector_tile_map.parquet', 'overlay_manifest.json']):
            parser.error('Overlay outputs already exist; use a new output directory.')
        vectors = load_vectors(args.vectors, geometry=True)
        tiles = raster_paths(args.rasters)
        metrics, mapping = overlay(vectors, tiles, args.sample_step)
        metrics.to_parquet(dest, index=False)
        mapping.to_parquet(args.output_dir / 'vector_tile_map.parquet', index=False)
        manifest = {'method': METHOD, 'sample_step_m': args.sample_step,
                    'vectors': str(args.vectors.resolve()), 'rasters': str(args.rasters.resolve()),
                    'vector_rows': len(vectors), 'tile_count': len(tiles), 'zones': ZONES,
                    'unassigned_vectors': int(metrics.dominant_cat.isna().sum()),
                    'warning': 'Legacy sampled length exceeds geometric length; coverage_pct is not true coverage. See docs/source_sink.md.'}
        (args.output_dir / 'overlay_manifest.json').write_text(json.dumps(manifest, indent=2))
        print(dest)


if __name__ == '__main__':
    main()
