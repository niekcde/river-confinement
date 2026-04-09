# river-confinement

Code for deriving river confinement metrics from SWORD river geometry and FABDEM elevation data.

## Status

The large refactor of the canonical pipeline is done.

For the active pipeline path:
- Steps `1.5` through `8` are now config-driven through `pipeline/paths.py`
- active code lives in [pipeline](/Users/6256481/Code/river-confinement/pipeline)
- legacy or superseded code lives in [old](/Users/6256481/Code/river-confinement/old)

What remains is mostly:
- performance work
- notebook cleanup
- reduction of unnecessary intermediates
- a few longer-term design changes tracked in [ISSUES.md](/Users/6256481/Code/river-confinement/ISSUES.md)

Resolved refactor work is tracked in [RESOLVED_ISSUES.md](/Users/6256481/Code/river-confinement/RESOLVED_ISSUES.md).

## Repository Layout

- [pipeline](/Users/6256481/Code/river-confinement/pipeline): active pipeline code and the final-results notebook
- [analysis](/Users/6256481/Code/river-confinement/analysis): additional analysis material
- [old](/Users/6256481/Code/river-confinement/old): legacy code kept for reference
- [config](/Users/6256481/Code/river-confinement/config): path config files
- [results](/Users/6256481/Code/river-confinement/results): default generated outputs
- [input_created](/Users/6256481/Code/river-confinement/input_created): derived DEM cache files such as the FABDEM VRT

## Environment

Create the environment with:

```bash
conda env create -f environment.yml
conda activate test-orthogonals
```

`pyarrow` is required for the canonical Step 2 bend-table output.

## Configuration

Create a local path config by copying [config/paths.example.json](/Users/6256481/Code/river-confinement/config/paths.example.json) to `config/paths.local.json`, or pass a config file explicitly on the CLI, or set `RIVER_CONFINEMENT_PATHS`.

Example:

```json
{
  "data_root": "./data",
  "results_root": "./results",
  "swot_vector_dir": "/Volumes/PhD/SWORD/v17/GPKG",
  "swot_nodes_dir": "/Volumes/PhD/SWORD/v17/GPKG",
  "fabdem_dir": "/Volumes/PhD/FABDEM",
  "input_created_root": "./input_created"
}
```

Important path behavior:
- `results_root` is the root for all generated pipeline outputs
- `input_created_root` is the root for derived DEM cache files
- the canonical pipeline does not require hard-coded `/Volumes/...` or `/scratch/...` paths in code anymore

## Canonical Pipeline

The canonical pipeline stages are:

1. Step `1.5`: build FABDEM VRT and bounds cache
2. Step `1`: segment SWORD reaches and build reference tables
3. Step `2`: build canonical bend-level Parquet tables
4. Step `4`: build the confinement-factor lookup table
5. Step `5`: compute per-file confinement outputs for a height factor
6. Step `6`: aggregate Step 5 outputs by height factor
7. Step `7`: spatially smooth the aggregated outputs
8. Step `8`: run clustering analysis on the smoothed outputs

Step `3` still exists, but only as a compatibility exporter from bend Parquet to the older CSV format. It is not required for the canonical path.

### Step 1.5

```bash
python -m pipeline.build_fabdem_index --config config/paths.local.json
```

Writes:
- `input_created/dem/FAB_dem_vrt.vrt`
- `input_created/dem/FAB_dem_bounds.gpkg`

### Step 1

```bash
python -m pipeline.segment_reaches af eu na sa as oc --workers 1 --config config/paths.local.json
```

Writes:
- `results/new_segments/vector/{file}_{group}_reach_new_segments.gpkg`
- `results/new_segments/node/{file}_{group}_node_new_segments.gpkg`
- `results/new_segments/vector_cont/{continent}_reaches.gpkg`
- `results/reference_tables/SWORD_stats.csv`
- `results/reference_tables/smoothingFactor.csv`
- `results/reference_tables/file_sorting.csv`

### Step 2

```bash
python -m pipeline.build_step2_results af 4 --config config/paths.local.json
```

For single-file debugging:

```bash
python -m pipeline.build_step2_results \
  --vector-file results/new_segments/vector/af_00_reach_new_segments.gpkg \
  --config config/paths.local.json
```

Writes:
- `results/orthogonals/{file}_50.gpkg`
- `results/bends/{file}_50.parquet`

The bend-level Parquet in `results/bends/` is the canonical downstream handoff.

### Step 4

```bash
python -m pipeline.build_step4_confinement_factor --config config/paths.local.json
```

Writes:
- `results/reference_tables/confinement_factor_50.csv`

### Step 5

```bash
python -m pipeline.build_step5_confinement_outputs af 4 --height-factor 2 --config config/paths.local.json
```

For single-file debugging:

```bash
python -m pipeline.build_step5_confinement_outputs \
  --input-file results/bends/af_00_50.parquet \
  --height-factor 2 \
  --config config/paths.local.json
```

Writes:
- `results/reach_averaged/{file}_50_{hf}.gpkg`
- `results/single_values/{file}_50_{hf}_conf.gpkg`
- `results/single_values/{file}_50_{hf}_conf.nc`

### Step 6

```bash
python -m pipeline.build_step6_aggregates --height-factor 2 --config config/paths.local.json
```

Writes:
- `results/single_values/global_50_02_conf.nc`
- `results/reach_averaged/{continent}_50_02_reachAveraged_conf.gpkg`

### Step 7

```bash
python -m pipeline.build_step7_spatial_smoothing \
  --height-factor 2 \
  --workers 1 \
  --continents af eu na sa as oc \
  --config config/paths.local.json
```

Writes:
- `results/single_smoothed/{continent}_50_{hf}_smoothed.nc`
- `results/single_smoothed/length_dict_{continent}.pkl`
- `results/single_smoothed/global_50_{hf}_smoothed.nc`

### Step 8

```bash
python -m pipeline.build_step8_confinement_clustering \
  --height-factors 2 3 4 \
  --clusters 2 3 \
  --sample-size 10 \
  --workers 1 \
  --config config/paths.local.json
```

Writes:
- `results/KmeanScores_confinement_{hf}.csv`
- `results/GMMSCORES_confinement_{hf}.csv`

### Combined Step 7-8 Workflow

```bash
python -m pipeline.build_step7_8_confinement_workflow \
  --height-factors 2 3 4 \
  --continents af eu na sa as oc \
  --smoothing-workers 1 \
  --clustering-workers 1 \
  --clusters 2 3 \
  --sample-size 10 \
  --config config/paths.local.json
```

## Helper Scripts

There are helper shell scripts in the repo root for smoke tests and larger runs:
- [run_oc00_smoke_test.sh](/Users/6256481/Code/river-confinement/run_oc00_smoke_test.sh)
- [run_continent_test.sh](/Users/6256481/Code/river-confinement/run_continent_test.sh)
- [run_world_test.sh](/Users/6256481/Code/river-confinement/run_world_test.sh)

These are developer convenience scripts, not the source of truth. The source of truth is the Python stage entrypoints listed above.

## Main Outputs

The most important canonical outputs are:
- `results/bends/*.parquet`: canonical Step 2 bend tables
- `results/reference_tables/confinement_factor_50.csv`: Step 4 lookup table
- `results/single_values/global_50_{hf}_conf.nc`: aggregated unsmoothed confinement outputs
- `results/single_smoothed/global_50_{hf}_smoothed.nc`: aggregated smoothed outputs used in later analysis

## Final Notebook

[pipeline/final_results.ipynb](/Users/6256481/Code/river-confinement/pipeline/final_results.ipynb) now uses the shared path setup and the canonical smoothed outputs, but it is still partly manual and still contains some exploratory analysis logic.

For notebook runs, either:
- set `RIVER_CONFINEMENT_PATHS` before launching Jupyter
- or edit the config selection cell in the notebook header

## What Is Still Open

The main remaining work is tracked in [ISSUES.md](/Users/6256481/Code/river-confinement/ISSUES.md). At the moment that is mainly:
- Step 2 performance
- Step 7 performance
- reducing large intermediate outputs
- deciding on a canonical saved cluster-assignment output
- cleaning up the final notebook further
- a future move toward one global Parquet input instead of current Step 1 source discovery
