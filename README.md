# river-morphology
Code for analysing SWORD data and creating a different river morpohlogy metrics

## Confinement pipeline audit

### Step 1: build segmented reach inputs and prerequisite tables

Source of truth for this step is the root-level code in `main.py`, `reach_definition.py`, and `support.py`.

The first pipeline step is the prerequisite-generation path inside `main.py`. The default CLI path later in `main.py` depends on outputs from this step, so `run_confinement_values_shell.py` is not the start of the pipeline.

Current Step 1 code path:
- `main.py` -> `run_create_new_reaches_main(...)` -> `create_new_reaches_main(...)`
- `create_new_reaches_main(...)` reads `input/SWOT_vector/{continent}*17*gpkg` and `input/SWOT_nodes/{continent}*17*gpkg`
- For each vector/node pair it calls `new_reach_definition(..., save=True)`

Step 1 outputs written by code:
- `results/new_segments/vector/{file}_{group}_reach_new_segments.gpkg`
- `results/new_segments/node/{file}_{group}_node_new_segments.gpkg`
- `results/SWORD_stats.csv`
- `results/smoothingFactor.csv`
- `results/file_sorting.csv`
- `results/new_segments/vector_cont/{continent}_reaches.gpkg`

Notes from the code audit:
- This first step is currently gated in `main.py` behind `create_new = False`
- There is no `.sh`, `.bash`, or `.zsh` pipeline entrypoint in this repository
- Multiple active scripts hard-code `directory = '/scratch/6256481/'`; that external path was not verified in this session
- `main.py` creates `results/new_segments/vector/` and `results/new_segments/node/`, but I do not see code here creating `results/new_segments/vector_cont/` before `create_continent_new_reach(...)` writes to it

### Step 2: run the first confinement-metric pass on segmented reaches

Source of truth for this step is the default CLI path in `main.py`.

Current Step 2 code path:
- `main.py` reads `continentInput = sys.argv[1]` and `number_of_processors = int(sys.argv[2])`
- It loads `results/file_sorting.csv`, filters rows for the requested continent, and builds `multiInput = [[filePath, 50], ...]`
- Under `if __name__ == '__main__':` it runs `Pool(number_of_processors).imap(main, multiInput)`
- `main(...)` opens the segmented reach/node files from Step 1, reads `results/smoothingFactor.csv`, opens `input_created/FAB_dem_vrt.vrt`, and computes bend/confinement metrics

Step 2 outputs written by code:
- `results/all/{continent}_{file_id}_50.csv`

Notes from the code audit:
- Step 2 depends on Step 1 outputs already existing, especially `results/new_segments/...` and `results/file_sorting.csv`
- The active code fixes `confFactor = 50` in `main.py`
- These `results/all/*.csv` files are consumed downstream by `open_to_single_apex.py` and `run_confinement_values_shell.py`

### Step 3: expand Step 2 reach files into single-bend value tables

Source of truth for this step is `create_apex_val_dataframe(...)` in `run_confinement_values.py`.

Current Step 3 code path:
- `open_to_single_apex.py` imports `create_apex_val_dataframe`
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
- This same Step 3 transformation is also duplicated inside the `createNewSingleVal == True` block in `run_confinement_values_shell.py`
- Step 3 depends on Step 2 outputs already existing in `results/all/`
- I do not see active code creating `results/single_values/` before Step 3 writes into it

### Step 4: build the global confinement-factor lookup table

Source of truth for this step is the `createNewFactor == True` block in `run_confinement_values_shell.py` together with `confinement_factor_single_values(...)` in `support.py`.

Current Step 4 code path:
- `run_confinement_values_shell.py` reads all `results/single_values/??_??_50.csv`
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
- I do not see a standalone Step 4 entrypoint outside `run_confinement_values_shell.py`
- This lookup table is consumed in `run_confinement_values.py` when later confinement values are computed

### Step 5: compute confinement outputs for each height factor

Source of truth for this step is the `for hf in heightFactor` loop in `run_confinement_values_shell.py` together with `calc_confinement_values(...)` and `ER_slope_margin_values(...)` in `run_confinement_values.py`.

Current Step 5 code path:
- `run_confinement_values_shell.py` loops over `heightFactor = [2, 0.5, 1, 1.5, 3, 4, 6, 8, 10, 15]`
- For each `hf`, it builds `multiInput = [[singleValueFile, hf], ...]`
- Under `if __name__ == '__main__':` it parallelizes `run(...)`
- `run(...)` reads one `results/single_values/{continent}_{file_id}_50.csv` file and calls `calc_confinement_values(...)`

Transformation performed by code:
- `calc_confinement_values(...)` opens `input_created/FAB_dem_vrt.vrt`
- When reading CSV input, it converts nested list columns such as `distOut`, `distInn`, `elevInn`, and `elevOut` back from strings
- `ER_slope_margin_values(...)` reads `results/confinement_factor.csv`, assigns the nearest `conFactor` by `bendWidths`, and computes bend-scale confinement outputs
- The code then derives left/right ER and slope values from `LROrthog`
- The output is saved once per input file and once per `heightFactor`

Step 5 outputs written by code:
- `results/reach_averaged/{continent}_{file_id}_50_{hfSave}.gpkg`
- `results/single_values/{continent}_{file_id}_50_{hfSave}_conf.gpkg`
- `results/single_values/{continent}_{file_id}_50_{hfSave}_conf.nc`

Notes from the code audit:
- Step 5 depends on Step 3 outputs, Step 4 output `results/confinement_factor.csv`, and `input_created/FAB_dem_vrt.vrt`
- The shell script immediately runs aggregation helpers after each `heightFactor`, but I am treating those as the next pipeline step
- I do not see active code in the current audited flow that creates `input_created/FAB_dem_vrt.vrt`

### Step 6: aggregate Step 5 outputs by height factor

Source of truth for this step is `concat_nc_conf_files(...)` and `concat_reachAveraged(...)` in `run_confinement_values.py`, called from the `for hf in heightFactor` loop in `run_confinement_values_shell.py`.

Current Step 6 code path:
- After each Step 5 batch for one `hf`, `run_confinement_values_shell.py` calls `concat_nc_conf_files(directory, crossFactor, hf)`
- It then calls `concat_reachAveraged(directory, crossFactor, hf)`

Transformation performed by code:
- `concat_nc_conf_files(...)` opens all per-file `*_conf.nc` outputs for a given `hf`, appends continent/file identifiers, drops several geometry-heavy variables, and concatenates them into one global dataset
- `concat_reachAveraged(...)` opens all per-file reach-averaged GeoPackages for a given `hf`, appends file identifiers, and concatenates them into one continent-level GeoPackage

Step 6 outputs written by code:
- `results/single_values/global_50_{hf}_conf.nc`
- `results/reach_averaged/{continent}_50_{hf}_reachAveraged_conf.gpkg`

Notes from the code audit:
- Step 6 depends on Step 5 outputs already existing for the requested `hf`
- This is the last clearly visible stage of the core confinement-output pipeline

### Step 7: spatially smooth the aggregated global confinement dataset

Source of truth for this step is `spatial_smoothing.py` together with `concat_nc_smooth_files(...)` in `support.py`.

Current Step 7 code path:
- `spatial_smoothing.py` hard-codes `cross = 50` and `hfList = [2]`
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
- This is downstream post-processing of Step 6 output, not part of the core confinement extraction itself
- In current code, only height factor `02` is processed here

### Step 8: run confinement clustering analysis on the smoothed dataset

Source of truth for this step is `clustering_confinement.py`.

Current Step 8 code path:
- `clustering_confinement.py` opens `results/single_smoothed/global_50_0{ch}_smoothed.nc`
- It loops over `ch in [2,3,4]`
- It scales selected smoothed slope variables, runs PCA, and then runs clustering-tuning code

Visible outputs written by code:
- `results/GMMSCORES_confinement_0{ch}.csv`

Notes from the code audit:
- This is downstream analysis, not a core confinement-generation step
- I do not see active KMeans output writes because the KMeans save line is commented out
- The script currently appears internally inconsistent with the upstream smoothing script, which only builds smoothed data for `hf = 02`
