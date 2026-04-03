# river-morphology
Code for analysing SWORD data and creating a different river morpohlogy metrics

## Planned Target Structure

This section is the working reference for the cleanup/refactor. It describes the intended stable project structure, not the claim that every part is already fully implemented.

Target layout:

```text
project/
  code/
    pipeline/
    analysis/
    old/
  config/
    paths.example.json
    paths.local.json
  input_created/
    dem/
      FAB_dem_vrt.vrt
      FAB_dem_bounds.gpkg
  results/
    reference_tables/
    new_segments/
      vector/
      node/
      vector_cont/
    orthogonals/
    profiles/
    all/
    single_values/
    reach_averaged/
    single_smoothed/
    logs/
    tmp/
  data/                     # optional
  README.md
  ISSUES.md
  environment.yml
```

Design rules for this target structure:
- `code/` holds repository code only
- `pipeline/` holds active pipeline and final-workflow code
- `analysis/` holds additional or currently unused analysis scripts and notebooks
- `old/` holds legacy or redundant code
- `results/` holds all generated outputs from the pipeline
- `data/` is optional and should not be mandatory for running the code if shared project data already exists elsewhere

Data-location principle:
- The code should support both a local project `data/` folder and an external shared data location
- The preferred long-term mechanism is a small local config file where the user sets a few root paths once
- The pipeline should then derive the rest of its paths from that config instead of hard-coded absolute paths

Planned path-config idea:
- `config/paths.example.json` documents the expected keys
- `config/paths.local.json` is user-specific and not intended as a shared hard-coded project path file
- Typical keys will include a results root plus shared-data roots such as SWORD vectors, SWORD nodes, and FABDEM

Current implementation status:
- The canonical Steps 1 through 8 now use the shared path-config layer in `pipeline/paths.py`
- The legacy Step 3 to Step 6 wrapper in `pipeline/run_confinement_values_shell.py` is now path-aware and delegates to the canonical stage entrypoints
- `pipeline/select_raster.py` and the active path-setup cells in `pipeline/final_results.ipynb` now also use the shared project paths
- `config/paths.local.json` is the intended user-specific override file and is gitignored
- Remaining path cleanup is now mostly limited to older non-canonical helper functions, compatibility arguments, and stale commented examples

Implementation note:
- The physical top-level organization may move toward `code/`
- The importable Python package should remain `pipeline`, not `code`, to avoid conflicts with Python's standard library `code` module

## Confinement pipeline audit

### Step 1: build segmented reach inputs and prerequisite tables

Source of truth for this step is the active pipeline code in `pipeline/segment_reaches.py`, `pipeline/paths.py`, `pipeline/reach_definition.py`, and `pipeline/support.py`.

The first pipeline step is now a dedicated segmentation entrypoint. The default CLI path later in `pipeline/main.py` depends on outputs from this step, so `pipeline/run_confinement_values_shell.py` is not the start of the pipeline.

Current Step 1 code path:
- Run `python -m pipeline.segment_reaches [continents...]`
- `pipeline/segment_reaches.py` -> `run_step1_segmentation(...)` -> `segment_continent(...)`
- `pipeline/paths.py` loads `config/paths.local.json` when present, otherwise falls back to local `data/` or legacy `input/` locations
- `segment_continent(...)` reads `{swot_vector_dir}/{continent}*17*gpkg` and `{swot_nodes_dir}/{continent}*17*gpkg`
- For each vector/node pair it calls `new_reach_definition(..., save=True)`
- `pipeline/segment_reaches.py` -> `build_step1_reference_tables(...)` rebuilds the prerequisite CSVs and continent-level segmented outputs

Step 1 outputs written by code:
- `results/new_segments/vector/{file}_{group}_reach_new_segments.gpkg`
- `results/new_segments/node/{file}_{group}_node_new_segments.gpkg`
- `results/reference_tables/SWORD_stats.csv`
- `results/reference_tables/smoothingFactor.csv`
- `results/reference_tables/file_sorting.csv`
- `results/new_segments/vector_cont/{continent}_reaches.gpkg`

Notes from the code audit:
- This first step now has its own stable entrypoint in `pipeline/segment_reaches.py`
- `pipeline/main.py` now starts at Step 2
- There is no `.sh`, `.bash`, or `.zsh` pipeline entrypoint in this repository
- Multiple active scripts hard-code `directory = '/scratch/6256481/'`; that external path was not verified in this session
- The current Step 1 entrypoint still assumes exactly one reach file and one node file per continent, because the downstream Step 2 filename logic expects outputs named like `{continent}_{group}_...`
- The Step 1 path/config layer is now separate from the Step 2+ hard-coded path handling, so the rest of the pipeline still needs the same cleanup pattern later

### Step 1.5: build the FABDEM VRT and bounds cache

Source of truth for this step is the active pipeline code in `pipeline/build_fabdem_index.py`, `pipeline/paths.py`, and `pipeline/dem.py`.

Current Step 1.5 code path:
- Run `python -m pipeline.build_fabdem_index`
- `pipeline/build_fabdem_index.py` loads the configured FABDEM source directory from `config/paths.local.json`
- It collects all FABDEM GeoTIFF tiles from `{fabdem_dir}`
- It builds the derived DEM cache files used later in the pipeline

Step 1.5 outputs written by code:
- `input_created/dem/FAB_dem_vrt.vrt`
- `input_created/dem/FAB_dem_bounds.gpkg`

Notes from the code audit:
- This is now an explicit intermediate pipeline step instead of an undeclared prerequisite
- The config is intended to store the source FABDEM directory, while the VRT and bounds paths are derived by the pipeline itself
- Step 2 and later active DEM-dependent code now read the centralized `input_created/dem/FAB_dem_vrt.vrt` path

### Step 2: run the first confinement-metric pass on segmented reaches

Source of truth for this step is the default CLI path in `pipeline/main.py`.

Current Step 2 code path:
- `pipeline/main.py` reads `continentInput = sys.argv[1]` and `number_of_processors = int(sys.argv[2])`
- It loads `results/reference_tables/file_sorting.csv`, filters rows for the requested continent, and builds `multiInput = [[filePath, 50], ...]`
- Under `if __name__ == '__main__':` it runs `Pool(number_of_processors).imap(main, multiInput)`
- `main(...)` opens the segmented reach/node files from Step 1, reads `results/reference_tables/smoothingFactor.csv`, and computes bend/confinement metrics plus orthogonal geometry
- It writes a geometry intermediate to `results/orthogonals/{continent}_{file_id}_50.gpkg`
- The geometry-only stage can now be run directly with `python -m pipeline.build_step2_orthogonals {continent} {processors}`
- For single-file debugging, the geometry-only stage can also be run with `python -m pipeline.build_step2_orthogonals --vector-file results/new_segments/vector/{file}.gpkg`
- It then runs the separate DEM-sampling stage in `pipeline/sample_step2_profiles.py`
- That stage reads the saved orthogonals, opens `input_created/dem/FAB_dem_vrt.vrt`, and writes sampled profiles to `results/profiles/{continent}_{file_id}_50.csv`
- `pipeline/main.py` finally assembles those sampled values back into the final Step 2 table

Step 2 outputs written by code:
- `results/orthogonals/{continent}_{file_id}_50.gpkg`
- `results/profiles/{continent}_{file_id}_50.csv`
- `results/all/{continent}_{file_id}_50.csv`

Notes from the code audit:
- Step 2 depends on Step 1 outputs already existing, especially `results/new_segments/...` and `results/reference_tables/file_sorting.csv`
- Step 2 also depends on Step 1.5 output `input_created/dem/FAB_dem_vrt.vrt`
- The active code fixes `confFactor = 50` in `pipeline/main.py`
- The current Step 2 implementation now pre-groups reach and node rows per `combined_reach_id` inside `pipeline/main.py` to reduce repeated dataframe filtering within a file
- The current Step 2 implementation now separates orthogonal-line geometry creation from DEM profile sampling inside `pipeline/get_orthogonals.py`
- `pipeline/main.py` now writes a real orthogonal intermediate and assembles the final Step 2 CSV from a separate sampled-profile intermediate
- The new geometry-only entrypoint in `pipeline/build_step2_orthogonals.py` allows Step 2 orthogonals to be tested independently from DEM sampling
- The DEM-sampling stage can now also be called directly with `python -m pipeline.sample_step2_profiles results/orthogonals/{file}.gpkg`
- When Step 2 is run with `processors = 1`, the current code now runs sequentially instead of using a one-worker multiprocessing pool
- These `results/all/*.csv` files are consumed downstream by `pipeline/open_to_single_apex.py` and `pipeline/run_confinement_values_shell.py`

### Step 3: expand Step 2 reach files into single-bend value tables

Source of truth for this step is `pipeline/build_step3_single_values.py` together with `create_apex_val_dataframe(...)` in `pipeline/run_confinement_values.py`.

Current Step 3 code path:
- Run `python -m pipeline.build_step3_single_values {continent} {processors}`
- For single-file debugging, run `python -m pipeline.build_step3_single_values --input-file results/all/{file}.csv`
- `pipeline/build_step3_single_values.py` discovers Step 2 files in `results/all/`, ensures `results/single_values/` exists, and parallelizes the Step 3 worker
- The worker calls `create_apex_val_dataframe(...)`
- `pipeline/open_to_single_apex.py` is now only a compatibility wrapper around the canonical Step 3 CLI
- The Step 3 block in `pipeline/run_confinement_values_shell.py` now also calls the same helper instead of duplicating its own file loop

Transformation performed by code:
- `create_apex_val_dataframe(...)` reads one Step 2 CSV
- It keeps only rows where `include_flag == '0'` and `calculated == '0'`
- It drops several reach-level columns and converts stored string/list fields back to Python objects
- It expands list-valued bend fields into one row per bend with `expand_dataframe(...)`

Step 3 outputs written by code:
- `results/single_values/{continent}_{file_id}_50.csv`

Notes from the code audit:
- Step 3 depends on Step 2 outputs already existing in `results/all/`
- The canonical Step 3 path is now config-driven through `pipeline/paths.py`
- `results/single_values/` is now created explicitly by the Step 3 helper before writes
- The Step 3 parser now has to normalize Step 2 nested profile columns that are stored as stringified numpy-style list values in `results/all/*.csv`

### Step 4: build the global confinement-factor lookup table

Source of truth for this step is `pipeline/build_step4_confinement_factor.py` together with `confinement_factor_single_values(...)` in `pipeline/support.py`.

Current Step 4 code path:
- Run `python -m pipeline.build_step4_confinement_factor`
- For a continent subset, run `python -m pipeline.build_step4_confinement_factor {continent}`
- For single-file debugging, run `python -m pipeline.build_step4_confinement_factor --input-file results/single_values/{file}.csv`
- `pipeline/build_step4_confinement_factor.py` discovers Step 3 files in `results/single_values/`
- It reads only the `bendWidths` column from those files
- It calls `confinement_factor_single_values(df, 'bendWidths', 50, 10)`
- It writes the result to `results/reference_tables/confinement_factor_50.csv`
- The Step 4 block in `pipeline/run_confinement_values_shell.py` now also calls the same helper instead of keeping its own concatenation logic

Transformation performed by code:
- `confinement_factor_single_values(...)` computes a linearly scaled factor from `bendWidths`
- It then inverts that value with `1 / cf`
- The saved table contains `bendWidths` and `conFactor`

Step 4 outputs written by code:
- `results/reference_tables/confinement_factor_50.csv`

Notes from the code audit:
- Step 4 depends on Step 3 outputs already existing in `results/single_values/`
- The canonical Step 4 path is now config-driven through `pipeline/paths.py`
- The new Step 4 helper now fails with a clear file error when no Step 3 inputs are available
- This lookup table is consumed in `pipeline/run_confinement_values.py` when later confinement values are computed

### Step 5: compute confinement outputs for each height factor

Source of truth for this step is `pipeline/build_step5_confinement_outputs.py` together with `calc_confinement_values(...)` and `ER_slope_margin_values(...)` in `pipeline/run_confinement_values.py`.

Current Step 5 code path:
- Run `python -m pipeline.build_step5_confinement_outputs {continent} {processors} --height-factor {hf}`
- For single-file debugging, run `python -m pipeline.build_step5_confinement_outputs --input-file results/single_values/{file}.csv --height-factor {hf}`
- `pipeline/build_step5_confinement_outputs.py` discovers Step 3 files in `results/single_values/`
- It ensures `results/reach_averaged/` and the Step 5 `results/single_values/*_conf.*` output paths exist
- The worker reads one `results/single_values/{continent}_{file_id}_50.csv` file and calls `calc_confinement_values(...)`
- The `for hf in heightFactor` loop in `pipeline/run_confinement_values_shell.py` now also calls the same helper instead of keeping its own per-file multiprocessing block

Transformation performed by code:
- `calc_confinement_values(...)` opens `input_created/dem/FAB_dem_vrt.vrt`
- When reading CSV input, it converts nested list columns such as `distOut`, `distInn`, `elevInn`, and `elevOut` back from strings
- `ER_slope_margin_values(...)` reads `results/reference_tables/confinement_factor_50.csv` through the shared path setup, assigns the nearest `conFactor` by `bendWidths`, and computes bend-scale confinement outputs
- The code then derives left/right ER and slope values from `LROrthog`
- The output is saved once per input file and once per `heightFactor`

Step 5 outputs written by code:
- `results/reach_averaged/{continent}_{file_id}_50_{hfSave}.gpkg`
- `results/single_values/{continent}_{file_id}_50_{hfSave}_conf.gpkg`
- `results/single_values/{continent}_{file_id}_50_{hfSave}_conf.nc`

Notes from the code audit:
- Step 5 depends on Step 3 outputs, Step 4 output `results/reference_tables/confinement_factor_50.csv`, and `input_created/dem/FAB_dem_vrt.vrt`
- The canonical Step 5 path is now config-driven through `pipeline/paths.py`
- `results/reach_averaged/` is now created explicitly by the Step 5 helper before writes
- The shell script immediately runs aggregation helpers after each `heightFactor`, but I am treating those as the next pipeline step

### Step 6: aggregate Step 5 outputs by height factor

Source of truth for this step is `pipeline/build_step6_aggregates.py` together with `concat_nc_conf_files(...)` and `concat_reachAveraged(...)` in `pipeline/run_confinement_values.py`.

Current Step 6 code path:
- Run `python -m pipeline.build_step6_aggregates --height-factor {hf}`
- `pipeline/build_step6_aggregates.py` calls `concat_nc_conf_files(...)` and `concat_reachAveraged(...)` for one `heightFactor`
- After each Step 5 batch for one `hf`, `pipeline/run_confinement_values_shell.py` now calls the same helper instead of calling the aggregation functions directly

Transformation performed by code:
- `concat_nc_conf_files(...)` opens all per-file `*_conf.nc` outputs for a given `hf`, appends continent/file identifiers, drops several geometry-heavy variables, and concatenates them into one global dataset
- `concat_reachAveraged(...)` opens all per-file reach-averaged GeoPackages for a given `hf`, appends file identifiers, and concatenates them into one continent-level GeoPackage

Step 6 outputs written by code:
- `results/single_values/global_50_{hf}_conf.nc`
- `results/reach_averaged/{continent}_50_{hf}_reachAveraged_conf.gpkg`

Notes from the code audit:
- Step 6 depends on Step 5 outputs already existing for the requested `hf`
- The canonical Step 6 path is now config-driven through `pipeline/paths.py`
- The Step 6 helper now formats the requested `hf` consistently and fails with a clear file error when no Step 5 inputs are available
- This is the last clearly visible stage of the unsmoothed confinement-output pipeline

### Step 7: spatially smooth the aggregated global confinement dataset

Source of truth for this step is `pipeline/build_step7_spatial_smoothing.py` together with `pipeline/spatial_smoothing.py` and `concat_nc_smooth_files(...)` in `pipeline/support.py`.

Current Step 7 code path:
- Run `python -m pipeline.build_step7_spatial_smoothing --height-factor {hf}`
- `pipeline/spatial_smoothing.py` reads `results/single_values/global_50_{hf}_conf.nc` through the shared path setup
- It derives additional normalized ER/slope variables from the global dataset
- It runs `run_bend_smoothing(...)` per continent, with optional continent filtering and configurable workers
- It writes one continent smoothed file per processed continent plus the corresponding `length_dict_{continent}.pkl`
- It then calls `concat_nc_smooth_files(...)` to build the global smoothed dataset

Transformation performed by code:
- The script builds a bend-neighbor graph within each continent/network
- It computes distance-weighted smoothed slope variables and smoothed standard deviations
- It writes one smoothed NetCDF per continent and one combined global smoothed NetCDF

Step 7 outputs written by code:
- `results/single_smoothed/{continent}_50_{hf}_smoothed.nc`
- `results/single_smoothed/length_dict_{continent}.pkl`
- `results/single_smoothed/global_50_{hf}_smoothed.nc`

Notes from the code audit:
- This is part of the main pipeline because the smoothed outputs feed the clustering workflow and the final-results workflow
- The canonical Step 7 path is now config-driven through `pipeline/paths.py`
- `results/single_smoothed/` is now created explicitly by the Step 7 path before writes

### Step 8: run confinement clustering analysis on the smoothed dataset

Source of truth for this step is `pipeline/build_step8_confinement_clustering.py` together with `pipeline/clustering_confinement.py`.

Current Step 8 code path:
- Run `python -m pipeline.build_step8_confinement_clustering --height-factors 2 3 4`
- `pipeline/clustering_confinement.py` reads `results/single_smoothed/global_50_{hf}_smoothed.nc` through the shared path setup
- It filters to rows with complete normalized slope and smoothed-slope variables
- It scales those variables, runs PCA, and then runs KMeans and GMM tuning over the requested cluster counts
- It writes one KMeans score table and one GMM score table per requested height factor

Visible outputs written by code:
- `results/KmeanScores_confinement_{hf}.csv`
- `results/GMMSCORES_confinement_{hf}.csv`

Notes from the code audit:
- This is part of the final analysis pipeline built on top of the smoothed confinement outputs
- The canonical Step 8 path is now config-driven through `pipeline/paths.py`
- The canonical Step 8 path now writes both KMeans and GMM tuning tables and no longer depends on the old undefined `comb` variable
- The clustering workflow still depends on Step 7 having already produced the requested smoothed height factors
- `pipeline/final_results.ipynb` imports clustering helpers and reads the smoothed `global_50_02_smoothed.nc` output, so the smoothing stage is not optional for the current final workflow
- The active path-setup cells in `pipeline/final_results.ipynb` now point at the shared project paths instead of the older `/Volumes/...` layout
