"""Reproducible companion to global_confinement_mollweide.ipynb.

Run with a Python environment containing geopandas, shapely >= 2, plotly,
numpy, pandas and pyproj. All outputs are local; the source is read-only.
"""
# %% Configuration
from pathlib import Path
import json
import sqlite3
import time
import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
import plotly.graph_objects as go
from matplotlib.colors import to_hex
import sys

SOURCE = Path('/Volumes/PhD/confinement/results/submission/v1/global_clusters_bend.gpkg')
LAND = Path.home() / 'Downloads/10m_physical.zip'
ROOT = Path.cwd() if (Path.cwd() / 'analysis').is_dir() else Path.cwd().parent
OUT = ROOT / 'results/global_confinement'
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / 'analysis'))
import importlib
import confinement_rendering
# Refresh the helper when this configuration is rerun in a notebook kernel.
importlib.reload(confinement_rendering)
from confinement_rendering import compact_map, land_background, OCEAN_COLOR, notebook_colormap, plotly_colorscale
MOLL = '+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs'
TOLERANCE = 250  # projected map metres: for display, not measurement
CAP = 50.0
COLORMAP = notebook_colormap(ROOT / 'analysis/global_confinement_mollweide.ipynb')
COLORS = [to_hex(COLORMAP(t)) for t in (0.0, 0.5, 1.0)]
PLOTLY_COLORSCALE = plotly_colorscale(COLORMAP)

# %% Load and audit
started = time.time()
gdf = gpd.read_file(SOURCE, layer='global_clusters_bend',
                    columns=['bendID', 'ER_inn', 'ER_out'])
assert gdf.crs is not None
# Strict arithmetic mean: one missing side => missing average.
ratio = (gdf.ER_inn + gdf.ER_out) / 2
finite = np.isfinite(ratio)
valid = finite & (ratio >= 0)
missing = ratio.isna()
invalid = ~finite & ~missing | (finite & (ratio < 0))
assert not invalid.any(), 'Review negative or infinite ratios before mapping.'
classes = np.where(ratio < 2, 0, np.where(ratio < 15, 1, 2))
stats = {
    'rows': len(gdf), 'source_bytes': SOURCE.stat().st_size,
    'missing_inner': int(gdf.ER_inn.isna().sum()),
    'missing_outer': int(gdf.ER_out.isna().sum()),
    'one_side_missing': int((gdf.ER_inn.isna() ^ gdf.ER_out.isna()).sum()),
    'missing_average': int(missing.sum()), 'above_50': int((ratio > CAP).sum()),
    'zeros': int((ratio == 0).sum()),
    'missing_or_empty_geometry': int((gdf.geometry.isna() | gdf.geometry.is_empty).sum()),
}
quantiles = ratio[valid].describe(percentiles=[.01, .1, .25, .5, .75, .9, .95, .99])
counts = pd.DataFrame({'class': ['0 ≤ r < 2', '2 ≤ r < 15', '15–50 (+ >50 capped, NaN)'],
                      'count': np.bincount(classes, minlength=3)})
counts['percent_all_bends'] = 100 * counts['count'] / len(gdf)
print(json.dumps(stats, indent=2), '\n', quantiles, '\n', counts, flush=True)
counts.to_csv(OUT / 'class_counts.csv', index=False)
quantiles.to_csv(OUT / 'ratio_quantiles.csv')

# %% Simplification benchmark
# Mollweide tolerances measure visual displacement in the projected map, not
# uniform ground distance. Keep all original attributes/geometry in the source.
projected = gdf.to_crs(MOLL)
geoms = projected.geometry.array
base_vertices = int(shapely.get_num_coordinates(geoms).sum())
base_wkb = int(sum(map(len, shapely.to_wkb(geoms))))
bench = []
for tolerance in [0, 100, 250, 500, 1000]:
    simplified = shapely.simplify(geoms, tolerance, preserve_topology=True) if tolerance else geoms
    vertices = int(shapely.get_num_coordinates(simplified).sum())
    wkb = int(sum(map(len, shapely.to_wkb(simplified))))
    bench.append(dict(tolerance_map_m=tolerance, vertices=vertices,
                      vertex_reduction_pct=100 * (1 - vertices / base_vertices),
                      geometry_wkb_MB=wkb / 1e6,
                      geometry_reduction_pct=100 * (1 - wkb / base_wkb)))
    if tolerance == TOLERANCE:
        display_geoms = simplified
benchmark = pd.DataFrame(bench)
benchmark.to_csv(OUT / 'simplification_benchmark.csv', index=False)
# Write paired representative samples with ALL attributes, so container and
# attribute overhead are included. This is an estimate, not a full-file claim.
sample = gpd.read_file(SOURCE, layer='global_clusters_bend',
                       where='fid % 100 = 0')
sample_simple = sample.to_crs(MOLL)
sample_simple.geometry = sample_simple.geometry.simplify(TOLERANCE, preserve_topology=True)
sample_simple = sample_simple.to_crs(sample.crs)
sample_sizes = {}
for name, frame in [('original', sample), ('simplified_250m', sample_simple)]:
    path = OUT / f'benchmark_sample_{name}.gpkg'
    if path.exists():
        path.unlink()
    frame.to_file(path, driver='GPKG', layer='bends')
    sample_sizes[name] = path.stat().st_size
stats['sample_rows'] = len(sample)
stats['sample_gpkg_bytes'] = sample_sizes
stats['sample_file_reduction_pct'] = 100 * (1 - sample_sizes['simplified_250m'] / sample_sizes['original'])
print(benchmark.to_string(index=False), '\nSample file sizes:', sample_sizes, flush=True)

# %% Map utilities: true projected linework with offline Natural Earth background
land = gpd.read_file(f'zip://{LAND}!ne_10m_land.shp').to_crs(MOLL)

def line_xy(geometries):
    """Flatten line parts with NaN separators, retaining every bend."""
    parts = shapely.get_parts(np.asarray(geometries, dtype=object))
    parts = parts[~shapely.is_empty(parts) & ~shapely.is_missing(parts)]
    if not len(parts):
        return [], []
    coords, ids = shapely.get_coordinates(parts, return_index=True)
    positions = np.arange(len(coords)) + ids
    xy = np.full((len(coords) + len(parts), 2), np.nan)
    xy[positions] = coords
    return xy[:, 0], xy[:, 1]

def base_map(title):
    fig = go.Figure()
    # Separate polygons prevent spurious bridges between islands. Natural Earth
    # land exteriors are sufficient here; inland lakes are not rendered.
    fig.add_trace(land_background(land))
    fig.update_layout(title=title, template='plotly_white', width=1400, height=780,
                      paper_bgcolor=OCEAN_COLOR, plot_bgcolor=OCEAN_COLOR,
                      margin=dict(l=15, r=100, t=85, b=65),
                      xaxis=dict(visible=False, range=[-18050000, 18050000]),
                      yaxis=dict(visible=False, range=[-9030000, 9030000],
                                 scaleanchor='x', scaleratio=1),
                      legend=dict(x=0.01, y=.12, bgcolor='rgba(255,255,255,.85)'),
                      annotations=[dict(x=0, y=-.065, xref='paper', yref='paper',
                         text='Mollweide • Natural Earth 1:10m land • each line is a bend • ratio = (inner + outer) / 2',
                         showarrow=False, font=dict(size=11))])
    return fig

def add_lines(fig, mask, color, name, showlegend=True):
    x, y = line_xy(display_geoms[np.asarray(mask)])
    fig.add_trace(go.Scattergl(x=x, y=y, mode='lines', name=name,
                              line=dict(color=color, width=1),
                              hoverinfo='name', showlegend=showlegend))

def save_map(fig, name):
    # Full vectors remain available separately. Default maps rasterize every
    # bend into a compact layer, avoiding million-line browser rendering.
    if name.startswith('mollweide_'):
        fig.write_html(OUT / f'{name}_vectors.html', include_plotlyjs=True)
        compact = compact_map(fig, OUT / f'{name}.png')
        fig.data = []
        fig.add_traces(compact.data)
        fig.layout = compact.layout
        fig.write_json(OUT / f'{name}_compact.json')
    fig.write_html(OUT / f'{name}.html', include_plotlyjs=True)
    print('Saved', OUT / f'{name}.html', flush=True)

# %% Map 1: requested classes
categorical = base_map('Average inner/outer bend confinement ratio — requested classes')
# Draw high/missing first, so low-valued confined bends remain visible on top.
labels = ['0–2', '2–15', '15–50 (also >50 and NaN)']
for cls in [2, 1, 0]:
    add_lines(categorical, classes == cls, COLORS[cls], labels[cls])
categorical.update_layout(legend_traceorder='reversed')
save_map(categorical, 'mollweide_classes')

# %% Map 2: logarithmic colour scale
# Plotly line traces accept one colour each. 96 narrow colour bands approximate
# a continuous logarithmic colour scale, with a colourbar in original units.
positive = valid & (ratio > 0)
LOG_MIN = min(1., float(ratio[positive].min()))

def continuous_map(transform, inverse, name, title):
    is_log = name == 'mollweide_log'
    low = LOG_MIN if is_log else 0.
    mapped = np.clip(ratio.to_numpy(), low, CAP)
    t0, t1 = float(transform(low)), float(transform(CAP))
    normalized = (transform(mapped) - t0) / (t1 - t0)
    bands = np.floor(np.nan_to_num(normalized, nan=0) * 96).clip(0, 95).astype(int)
    fig = base_map(title)
    add_lines(fig, missing, COLORS[-1], 'NaN (same colour as upper class)')
    if is_log and (ratio == 0).any():
        add_lines(fig, ratio == 0, '#888888', 'Zero (undefined on log scale)')
    for band in reversed(range(96)):
        mask = (positive if is_log else valid) & (bands == band)
        if mask.any():
            color = to_hex(COLORMAP((band + .5) / 96))
            lo, hi = inverse(t0 + np.array([band, band + 1]) / 96 * (t1 - t0))
            add_lines(fig, mask, color, f'{lo:.3g}–{hi:.3g}', False)
    ticks = [v for v in [low, 1, 2, 5, 10, 15, 25, 50] if low <= v <= CAP]
    ticks = sorted(set(ticks))
    fig.add_trace(go.Scattergl(x=[None, None], y=[None, None], mode='markers',
        marker=dict(color=[t0, t1], cmin=t0, cmax=t1, colorscale=PLOTLY_COLORSCALE,
                    showscale=True, colorbar=dict(title='Mean ratio',
                    tickvals=[float(transform(v)) for v in ticks],
                    ticktext=[f'{v:.3g}' if v < CAP else '50+' for v in ticks])),
        showlegend=False, hoverinfo='skip'))
    save_map(fig, name)
    return fig

logarithmic = continuous_map(np.log10, lambda x: 10 ** x, 'mollweide_log',
                             'Average inner/outer bend confinement ratio — logarithmic colours')
log1p_map = continuous_map(np.log1p, np.expm1, 'mollweide_log1p',
                           'Average inner/outer bend confinement ratio — log(1 + ratio) colours')

# %% Distribution and transformation comparison
# Compare where the two class boundaries fall when a transformed 0–50 scale
# is divided into equal thirds. Agreement measures class membership, not a
# perceptual validation of the continuous colour maps; NaNs excluded.
comparisons = []
transforms = {
    'Linear 0–50': np.array([CAP / 3, 2 * CAP / 3]),
    'Log10 (positive; floor=min(1,min positive))': 10 ** (np.log10(LOG_MIN) + np.array([1, 2]) / 3 * np.log10(CAP / LOG_MIN)),
    'Log1p 0–50': np.expm1(np.array([1, 2]) / 3 * np.log1p(CAP)),
    'Square root 0–50': CAP * (np.array([1, 2]) / 3) ** 2,
    'Cube root 0–50': CAP * (np.array([1, 2]) / 3) ** 3,
    'Empirical terciles': np.quantile(ratio[valid], [1 / 3, 2 / 3]),
    'Piecewise linear anchors 0,2,15,50': np.array([2., 15.]),
}
for name, edges in transforms.items():
    mask = positive if name.startswith('Log10') else valid
    inferred = np.digitize(ratio[mask], edges)
    comparisons.append(dict(scale=name, lower_break=edges[0], upper_break=edges[1],
                            agreement_percent=100 * np.mean(inferred == classes[mask])))
comparison = pd.DataFrame(comparisons)
comparison.to_csv(OUT / 'legend_comparison.csv', index=False)
hist_count, hist_edges = np.histogram(np.log10(ratio[positive]), bins=100)
distribution = go.Figure(go.Bar(x=(hist_edges[:-1] + hist_edges[1:]) / 2,
                               y=hist_count, width=np.diff(hist_edges),
                               marker_color=COLORS[1], name='Bends'))
for edge in [2, 15, 50]:
    distribution.add_vline(x=np.log10(edge), line_dash='dash', annotation_text=str(edge))
distribution.update_layout(template='plotly_white', title='Distribution of positive mean ratios (missing and zero excluded)',
                            xaxis_title='log10(mean confinement ratio)', yaxis_title='Number of bends')
save_map(distribution, 'ratio_distribution')
stats['log_color_floor'] = LOG_MIN
(OUT / 'audit.json').write_text(json.dumps(stats, indent=2))
report = f'''# Global bend confinement: inspection

Source: `{SOURCE}`. Arithmetic mean of ER_inn and ER_out; either missing side makes the mean missing.
Counts are per bend, not length weighted. Class intervals are [0,2), [2,15), and [15,50],
with >50 capped and all missing means displayed in the final class. Missing is a display
convention, not evidence of a high ratio. See audit.json for missingness and out-of-range counts.

```
{counts.to_string(index=False)}

{quantiles.to_string()}

{comparison.to_string(index=False)}

{benchmark.to_string(index=False)}
```

The sample GeoPackages retain all attributes ({len(sample):,} systematically selected bends).
Original sample: {sample_sizes['original']/1e6:.2f} MB; simplified sample:
{sample_sizes['simplified_250m']/1e6:.2f} MB; reduction {stats['sample_file_reduction_pct']:.1f}%.
Full geometry reductions above are measured on all bends; sample file reductions are estimates
of container savings, not a measured full-file reduction. Rewriting a file is needed to reclaim
space; editing geometries in place may leave free SQLite pages. Keep the original for analysis.
The map uses {TOLERANCE} projected metres of simplification, a display-space tolerance rather
than a uniform ground-distance guarantee. Natural Earth is already generalized.

For an exactly matching outlook, retain the discrete thresholds 2 and 15, even if the legend
is positioned on a logarithmic axis. A continuous piecewise scale anchored at 0,2,15,50 preserves
both break locations at one-third and two-thirds but changes colours within each class.
The table quantifies how closely ordinary transforms reproduce the requested membership.
Cube root is the recommended simple smooth alternative: equally spaced colour thirds on
a 0–50 cube-root scale fall at 1.85185 and 14.81481, very close to 2 and 15. Keep legend tick
labels in original ratio units (0, 2, 15, 50). This reproduces class areas closely, although
a continuous colour ramp necessarily introduces within-class variation.
Log1p handles zero and is another smooth alternative; the log10 map excludes zero from its
numeric colour scale. All continuous maps use 96 colour bands because Plotly lines have a
single colour per trace. Missing averages stay in a separate labelled trace with the high-class
colour, and values above 50 saturate the scale. Class agreement excludes missing averages.

Maps are self-contained Plotly HTML with locally projected Mollweide geometry and Natural Earth
land (no online basemap). Default maps rasterize ALL bend strokes into a 4320-pixel-wide layer;
there is no bend sampling. Land and legends remain vectors. Deep zoom magnifies raster pixels;
use the larger *_vectors.html files when geometry-level zoom is needed. Full-data PNGs are
embedded in the notebook, with compact Plotly figures displayed directly (no iframe URLs).
Inland lake holes are omitted from the background.
'''
(OUT / 'inspection.md').write_text(report)
print(comparison.to_string(index=False), '\nFinished in', round(time.time() - started), 'seconds', flush=True)
