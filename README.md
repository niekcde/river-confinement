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
- The initial path-config layer is now in place for Step 1 via `pipeline/paths.py`
- `config/paths.local.json` is the intended user-specific override file and is gitignored
- Later pipeline steps still need to be migrated onto the same config-driven path system

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
- `main(...)` opens the segmented reach/node files from Step 1, reads `results/reference_tables/smoothingFactor.csv`, opens `input_created/dem/FAB_dem_vrt.vrt`, and computes bend/confinement metrics

Step 2 outputs written by code:
- `results/all/{continent}_{file_id}_50.csv`

Notes from the code audit:
- Step 2 depends on Step 1 outputs already existing, especially `results/new_segments/...` and `results/reference_tables/file_sorting.csv`
- Step 2 also depends on Step 1.5 output `input_created/dem/FAB_dem_vrt.vrt`
- The active code fixes `confFactor = 50` in `pipeline/main.py`
- These `results/all/*.csv` files are consumed downstream by `pipeline/open_to_single_apex.py` and `pipeline/run_confinement_values_shell.py`

### Step 3: expand Step 2 reach files into single-bend value tables

Source of truth for this step is `create_apex_val_dataframe(...)` in `pipeline/run_confinement_values.py`.

Current Step 3 code path:
- `pipeline/open_to_single_apex.py` imports `create_apex_val_dataframe`
- It reads `results/all/??_??_50.csv`
- It removes existing `results/single_values/??_??_50.csv`
- Under `if __name__ == '__main__':` it runs `Pool(10).imap(create_apex_val_dataframe, files)`

Transformation performed by code:
- `create_apex_val_dataframe(...)` reads one Step 2 CSV
- It keeps only rows where `include_flag == '0'` and `calculated == '0'`
- It drops several reach-level columns and converts stored string/list fields back to Python objects
- It expands list-valued bend fields into one row per bend with `expand_dataframe(...)`

Step 3 outputs written by code:
- `results/single_values/{continent}_{file_id}_50.csv`

Notes from the code audit:
- This same Step 3 transformation is also duplicated inside the `createNewSingleVal == True` block in `pipeline/run_confinement_values_shell.py`
- Step 3 depends on Step 2 outputs already existing in `results/all/`
- I do not see active code creating `results/single_values/` before Step 3 writes into it

### Step 4: build the global confinement-factor lookup table

Source of truth for this step is the `createNewFactor == True` block in `pipeline/run_confinement_values_shell.py` together with `confinement_factor_single_values(...)` in `pipeline/support.py`.

Current Step 4 code path:
- `pipeline/run_confinement_values_shell.py` reads all `results/single_values/??_??_50.csv`
- It concatenates those files into one dataframe
- It calls `confinement_factor_single_values(dfA, 'bendWidths', 50, 10)`
- It writes the result to `results/confinement_factor.csv`

Transformation performed by code:
- `confinement_factor_single_values(...)` computes a linearly scaled factor from `bendWidths`
- It then inverts that value with `1 / cf`
- The saved table contains `bendWidths` and `conFactor`

Step 4 outputs written by code:
- `results/confinement_factor.csv`

Notes from the code audit:
- Step 4 depends on Step 3 outputs already existing in `results/single_values/`
- I do not see a standalone Step 4 entrypoint outside `pipeline/run_confinement_values_shell.py`
- This lookup table is consumed in `pipeline/run_confinement_values.py` when later confinement values are computed

### Step 5: compute confinement outputs for each height factor

Source of truth for this step is the `for hf in heightFactor` loop in `pipeline/run_confinement_values_shell.py` together with `calc_confinement_values(...)` and `ER_slope_margin_values(...)` in `pipeline/run_confinement_values.py`.

Current Step 5 code path:
- `pipeline/run_confinement_values_shell.py` loops over `heightFactor = [2, 0.5, 1, 1.5, 3, 4, 6, 8, 10, 15]`
- For each `hf`, it builds `multiInput = [[singleValueFile, hf], ...]`
- Under `if __name__ == '__main__':` it parallelizes `run(...)`
- `run(...)` reads one `results/single_values/{continent}_{file_id}_50.csv` file and calls `calc_confinement_values(...)`

Transformation performed by code:
- `calc_confinement_values(...)` opens `input_created/dem/FAB_dem_vrt.vrt`
- When reading CSV input, it converts nested list columns such as `distOut`, `distInn`, `elevInn`, and `elevOut` back from strings
- `ER_slope_margin_values(...)` reads `results/confinement_factor.csv`, assigns the nearest `conFactor` by `bendWidths`, and computes bend-scale confinement outputs
- The code then derives left/right ER and slope values from `LROrthog`
- The output is saved once per input file and once per `heightFactor`

Step 5 outputs written by code:
- `results/reach_averaged/{continent}_{file_id}_50_{hfSave}.gpkg`
- `results/single_values/{continent}_{file_id}_50_{hfSave}_conf.gpkg`
- `results/single_values/{continent}_{file_id}_50_{hfSave}_conf.nc`

Notes from the code audit:
- Step 5 depends on Step 3 outputs, Step 4 output `results/confinement_factor.csv`, and `input_created/dem/FAB_dem_vrt.vrt`
- The shell script immediately runs aggregation helpers after each `heightFactor`, but I am treating those as the next pipeline step

### Step 6: aggregate Step 5 outputs by height factor

Source of truth for this step is `concat_nc_conf_files(...)` and `concat_reachAveraged(...)` in `pipeline/run_confinement_values.py`, called from the `for hf in heightFactor` loop in `pipeline/run_confinement_values_shell.py`.

Current Step 6 code path:
- After each Step 5 batch for one `hf`, `pipeline/run_confinement_values_shell.py` calls `concat_nc_conf_files(directory, crossFactor, hf)`
- It then calls `concat_reachAveraged(directory, crossFactor, hf)`

Transformation performed by code:
- `concat_nc_conf_files(...)` opens all per-file `*_conf.nc` outputs for a given `hf`, appends continent/file identifiers, drops several geometry-heavy variables, and concatenates them into one global dataset
- `concat_reachAveraged(...)` opens all per-file reach-averaged GeoPackages for a given `hf`, appends file identifiers, and concatenates them into one continent-level GeoPackage

Step 6 outputs written by code:
- `results/single_values/global_50_{hf}_conf.nc`
- `results/reach_averaged/{continent}_50_{hf}_reachAveraged_conf.gpkg`

Notes from the code audit:
- Step 6 depends on Step 5 outputs already existing for the requested `hf`
- This is the last clearly visible stage of the unsmoothed confinement-output pipeline

### Step 7: spatially smooth the aggregated global confinement dataset

Source of truth for this step is `pipeline/spatial_smoothing.py` together with `concat_nc_smooth_files(...)` in `pipeline/support.py`.

Current Step 7 code path:
- `pipeline/spatial_smoothing.py` hard-codes `cross = 50` and `hfList = [2]`
- It reads `results/single_values/global_50_02_conf.nc`
- It derives additional normalized ER/slope variables from the global dataset
- It runs `run_bend_smoothing(...)` per continent
- It calls `concat_nc_smooth_files(directory, cross, hf)`

Transformation performed by code:
- The script builds a bend-neighbor graph within each continent/network
- It computes distance-weighted smoothed slope variables and smoothed standard deviations
- It writes one smoothed NetCDF per continent and one combined global smoothed NetCDF

Step 7 outputs written by code:
- `results/single_smoothed/{continent}_50_02_smoothed.nc`
- `results/single_smoothed/length_dict_{continent}.pkl`
- `results/single_smoothed/global_50_02_smoothed.nc`

Notes from the code audit:
- This is part of the main pipeline because the smoothed outputs feed the clustering workflow and the final-results workflow
- In current code, only height factor `02` is processed here

### Step 8: run confinement clustering analysis on the smoothed dataset

Source of truth for this step is `pipeline/clustering_confinement.py`.

Current Step 8 code path:
- `pipeline/clustering_confinement.py` opens `results/single_smoothed/global_50_0{ch}_smoothed.nc`
- It loops over `ch in [2,3,4]`
- It scales selected smoothed slope variables, runs PCA, and then runs clustering-tuning code

Visible outputs written by code:
- `results/GMMSCORES_confinement_0{ch}.csv`

Notes from the code audit:
- This is part of the final analysis pipeline built on top of the smoothed confinement outputs
- I do not see active KMeans output writes because the KMeans save line is commented out
- The script currently appears internally inconsistent with the upstream smoothing script, which only builds smoothed data for `hf = 02`
- `pipeline/final_results.ipynb` imports clustering helpers and reads the smoothed `global_50_02_smoothed.nc` output, so the smoothing stage is not optional for the current final workflow
