# Source–bypass–sink comparison

This analysis compares saved confinement classes with Martin & Lamb (2025,
Geology) routing domains. The dataset is a **categorical raster**, not a polygon
classification: 0 Ocean, 1 Sink, 2 Bypass, 3 Source, 4 No data.
The downloaded package's `QuickStart Instructions.txt` requests citation of both
the paper and [the dataset](https://figshare.com/s/660b1728e0ec2ff803c6).

The workflow lives in `analysis/source_sink.py` and
`analysis/source_sink_summary.py`. It is an optional analysis, not an additional
Step 1–8 pipeline stage. Existing environment.yml dependencies are sufficient.
Run commands from the repository root with the project Python environment.
No command runs at import time; summarizing never invokes raster sampling.

## Reuse the existing manuscript assignments

Place these **original paired files** in `data/source_sink/` (or pass their
locations explicitly):

- `vectors.parquet`: bend identities, GMM classification, elevations, lateral
  slopes, geometry, and original `vector_id`.
- `vector_raster_metrics.parquet`: original per-bend routing metrics.

These were found in
`RivWidth/Paper/newAnalysis/data/` in the author's OneDrive archive. Both have
1,490,573 rows. The ZIPs alone contain raster inputs, not bend assignments.
The original source code was `RivWidth/Paper/newAnalysis/source_sink.py` and
`source_sink.ipynb`; this integration refactors their sampling, category mapping,
and final plotting logic without importing the notebook's exploratory attempts.

```bash
python -m analysis.source_sink_summary \
  --vectors data/source_sink/vectors.parquet \
  --metrics data/source_sink/vector_raster_metrics.parquet \
  --output-dir results/source_sink --plot
```

Outputs include:

- `bend_assignments_topography.parquet`: classified bends with routing metrics.
- `group_summary.csv`: n, valid n, missing n, median, Q1, Q3, IQR width, and an
  n<300 flag, for all assignments and the zero-transition sensitivity cohort.
- `mann_whitney.csv`: two-sided asymptotic tests and signed rank-biserial effects.
- `regions.csv`: continent counts for Unconfined Source bends.
- `percentage_audit.csv`: Unconfined count and sampled-length shares.
- `figure_percentages.csv` and optional `source_sink.png`: all four confinement
  classes, with explicit sampled-length labeling.
- `audit.json`: input provenance, excluded classifications, and metric checks.

The join requires complete one-to-one `vector_id` agreement and unique
`(source_file, bendID)` identities. Never rebuild or reorder IDs to join an old
metrics table. Equal row counts alone do not establish correct pairing.
The historical GMM mapping is 0 Unconfined, 2 Confined, 5 Sym partial, and
1/3/4/6 Asym partial. This mapping is specific to the manuscript model; cluster
labels from a newly fitted model must be verified before using it here.

Elevation is `cp_height` (m). Lateral slope is the mean of `slope_inn` and
`slope_out` (m/m), requiring both; individual-side summaries are also saved.
Negative valid values are retained; -9999 and 99999 are treated as missing.
`cm_height` is an imposed confinement threshold, not an independent relief metric.
Lateral slopes contribute to class assignment, so their separation is not
independent topographic validation. Tests use disjoint populations (Source versus
remaining Unconfined, not Source versus the overlapping global population).
They do not correct for spatial dependence. Small/empty groups are flagged;
empty-group tests are not performed.

## Run a new overlay explicitly

### September 18, 2026 HF=2 rerun

The new fit has a different GMM label order. Its user-confirmed mapping is in
`config/source_sink_hf2_20260918.json`: Unconfined 1, Sym partial 3,
Asym partial 0/2/4/5, Confined 6. Always pass this mapping for the new fit.
The default mapping remains historical for reproducibility of the old tables.

The new smoothed file was moved to `results/single_smoothed/`; the old path under
`smoothing_revision/hf2_local_n3_alpha075/` is no longer its location.
The new preparation command joins the saved GMM assignments to the NetCDF on
all five saved identifiers, then joins all 43 geometry and bend-profile tables
on `(file, reach_id, combined_reach_id, apex, bendDistOut)`. It requires complete,
unique matches, records SHA-256 fingerprints of the inputs, and preserves the
global NetCDF index as `vector_id`. It does not use old routing assignments.

```bash
python -m analysis.source_sink_refresh metrics \
  --smoothed results/single_smoothed/global_50_02_smoothed.nc \
  --assignments results/smoothing_revision/hf2_gmm7/gmm7_hf2_assignments.parquet \
  --output-dir results/source_sink_hf2_20260918

python -m analysis.source_sink_refresh geometry \
  --metrics results/source_sink_hf2_20260918/new_metrics.parquet \
  --geometry-dir results/single_values --bends-dir results/bends \
  --output-dir results/source_sink_hf2_20260918

python -m analysis.source_sink overlay \
  --vectors results/source_sink_hf2_20260918/vectors.parquet \
  --rasters data/source_sink/mask_strat_241022.zip \
  --output-dir results/source_sink_hf2_20260918

python -m analysis.source_sink_summary \
  --vectors results/source_sink_hf2_20260918/vectors.parquet \
  --metrics results/source_sink_hf2_20260918/vector_raster_metrics.parquet \
  --class-map config/source_sink_hf2_20260918.json \
  --output-dir results/source_sink_hf2_20260918/summary --plot
```

Use a fresh output directory for each rerun. In the current machine's existing
environments, `geo` supports NetCDF preparation and summary, while `geo2`
supports geometry preparation and raster sampling. The declared project
environment contains both sets of dependencies.

When supplied, `bendHeight` and `lineSlope` are also summarized/tested. The
former is mean bend DEM elevation; `abs_lineSlope` is the magnitude of the
longitudinal fitted elevation gradient (m/m), **shared by bends in the sampled
combined reach**. Thus it provides a different topographic measure from the
lateral confinement slope, but bend-level observations are not independent.
Group summaries report each metric's valid and missing counts separately.

Raster sampling now caches each tile's band and vectorizes point interpolation.
Nearest-cell values, outside-tile/nodata behavior and endpoint weighting are
unchanged; edge tests and ten original bends reproduced the saved results.

Preparation splits geographic lines with longitude jumps exceeding 180 degrees
at the dateline before projection. The new HF2 run exposed eight such bends:
without splitting, projected lengths were about 24,000 km. Their source GPKGs
were left unchanged; only the overlay input geometry was corrected. The run's
`dateline_correction.csv` and `dateline_metrics_before.parquet` preserve the
before/after evidence. This correction is now applied by both preparation
commands so future reruns do not reproduce the wraparound artifact.

### Other new runs

This can be expensive globally. It is unnecessary for the existing manuscript
statistics. The following commands are separate, explicit operations:

```bash
python -m analysis.source_sink prepare \
  --vector-dir /path/to/clusters \
  --output results/source_sink_new/vectors.parquet

python -m analysis.source_sink overlay \
  --vectors results/source_sink_new/vectors.parquet \
  --rasters data/source_sink/mask_strat_241022.zip \
  --output-dir results/source_sink_new

python -m analysis.source_sink_summary \
  --vectors results/source_sink_new/vectors.parquet \
  --metrics results/source_sink_new/vector_raster_metrics.parquet \
  --output-dir results/source_sink_new/summary --plot
```

`prepare` reads sorted `??_??_bend_cluster.gpkg` files/layers, projects to
EPSG:8857 and assigns IDs for that new run. It refuses to overwrite a vector
table. `overlay` also accepts historical `vectors.parquet` (plain WKB explicitly
in EPSG:8857), or the GeoParquet produced by `prepare`. It reads mask TIFFs
directly inside the supplied ZIP; an extracted directory is also accepted.
It writes per-bend metrics, the tile mapping, and an overlay manifest, and refuses
to overwrite existing overlay outputs. Spatial indexing limits each tile's work
to intersecting vectors. Processing is sequential and deterministic; no DuckDB
extension downloads or global multiprocessing jobs are triggered.

## Historical numerical conventions

The default sample step is 50 m in EPSG:8857. For each clipped line part, the
original script sampled `ceil(length/50)+1` points including both endpoints and
gave **every point** a weight of `length/(n-1)`. Thus `length_m` is a historical
sample weight, generally greater than geometric line length, not exact intersected
length. Category percentages and entropy use these weights.

The historical `coverage_pct` expression effectively yields 1 wherever there are
valid samples. It is retained for reproducibility and must not be interpreted as
true coverage. Raster nodata is skipped; classes 0 and 4 are retained as explicit
categories. Transitions skip nodata gaps and are counted within clipped line
parts, not across tile/part boundaries. Dominant-category ties use the first
encountered category. Tile processing is now sorted, so cross-tile ties can differ
from the old asynchronous completion order. No-sample bends are retained with
missing dominant category, zero length, and zero coverage.

The figure selects `transitions == 0`, valid confinement class, and dominant
zone in Sink/Bypass/Source. In the original paired data, Unconfined has
**n=800,489**, of which **n=271,278** are Source: **27.9694598664% of sampled
length** versus **33.8890% of bends** (see CSV for full precision).
This is the origin of the reported 28%; it is not the unfiltered bend-count share.

For historical results, reuse the saved paired tables. This implementation
preserves the known sampling conventions rather than silently replacing them
with an exact line/raster intersection algorithm. Any future correction should
use a separate method version and regenerate its own outputs.

## Checks

```bash
python -m unittest discover -s tests -p 'test_source_sink.py'
```

Tests cover a known three-zone raster, endpoint weights, nodata transitions,
unmatched bends, non-contiguous IDs, invalid joins, and GPKG-to-GeoParquet
preparation. The integrated summary was compared numerically against the
original saved-table summary. A ten-bend real-data check against the downloaded
mask ZIP reproduced every saved metric (rtol=1e-10, atol=1e-8), without a global
overlay rerun. Raster tests ran in the local `geo2` environment; the local `geo`
environment currently has a missing GDAL/netCDF shared library, but can run
the saved-table summary because raster libraries are imported only on demand.
Raw ZIPs, paired Parquet inputs and generated outputs stay in the already ignored
`data/` and `results/` folders; the analysis code, documentation and tests belong
in Git.
