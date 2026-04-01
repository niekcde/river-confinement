# Known Issues

## 1. Hard-coded base paths

Current active scripts hard-code `directory = '/scratch/6256481/'` instead of taking a repo-relative path, config value, or CLI argument.

Impact:
- The pipeline is not portable as checked into this repository
- Running code from this session cannot assume that `/scratch/6256481/` exists

Examples in current active code:
- `pipeline/main.py`
- `pipeline/open_to_single_apex.py`
- `pipeline/run_confinement_values_shell.py`
- `pipeline/run_confinement_values.py` default argument
- `pipeline/select_raster.py`
- `pipeline/spatial_smoothing.py`
- `pipeline/clustering_confinement.py`

## 2. Step 1 entrypoint separation

Status update:
- Initial fix implemented: Step 1 now has a dedicated entrypoint in `pipeline/segment_reaches.py`
- Step 1 paths are now loaded through `pipeline/paths.py` and `config/paths.local.json`

Remaining follow-up:
- Move the later stages onto the same config-driven path system so Step 1 and Step 2+ use one consistent project layout

Impact that remains:
- Step 1 is now cleanly separated, but the rest of the pipeline still mixes the old hard-coded layout with the new Step 1 path configuration

## 15. Step 1 still assumes one reach file and one node file per continent

The new `pipeline/segment_reaches.py` entrypoint currently fails if more than one `*17*.gpkg` reach file or node file is found for a continent.

Why this exists:
- The downstream active pipeline still assumes Step 1 outputs are named like `{continent}_{group}_...`
- That naming scheme comes from `new_reach_definition(...)`, where `group` is the segmented network-group id used later by `pipeline/main.py`

Impact:
- Step 1 is now stable and explicit, but it is not yet generalized to multi-file-per-continent source layouts
- Supporting multiple source files per continent will require aligning Step 2 and the later file-sorting logic, not just Step 1

Follow-up direction:
- Keep the current one-file-per-continent assumption for now
- Revisit this only when the later stages are refactored onto the same path/file naming system

## 3. Step 3 has duplicated entrypoints

The Step 3 single-bend expansion is implemented once in `create_apex_val_dataframe(...)`, but it is invoked by both `pipeline/open_to_single_apex.py` and the `createNewSingleVal == True` block in `pipeline/run_confinement_values_shell.py`.

Impact:
- There are two current ways to run the same pipeline stage
- Pipeline documentation and execution order are harder to keep consistent

Follow-up direction:
- Keep one canonical Step 3 entrypoint and make the other call it or remove the duplication later

## 4. Downstream output directories are not created in active setup code

Active downstream code writes to `results/single_values/` and `results/reach_averaged/`, but the visible setup code in `pipeline/main.py` only creates:
- `results/`
- `results/new_segments/`
- `results/new_segments/node/`
- `results/new_segments/vector/`
- `results/all/`
- `results/centerline/`
- `results/cycles/`

Impact:
- Step 3 and later stages may fail on a clean filesystem if those output directories do not already exist

Follow-up direction:
- Add explicit directory creation for downstream outputs in the active pipeline setup path

## 5. Step 4 is only embedded in `pipeline/run_confinement_values_shell.py`

The confinement-factor lookup table is produced by the `createNewFactor == True` block in `pipeline/run_confinement_values_shell.py`. I do not see a separate dedicated entrypoint for this stage.

Impact:
- Step 4 is harder to run independently and document cleanly
- The same script mixes Step 3, Step 4, and later confinement-value stages

Follow-up direction:
- Split the lookup-table build into its own explicit stage entrypoint later

## 6. Step 4 is not robust to an empty Step 3 result set

In `pipeline/run_confinement_values_shell.py`, `dfA` is only assigned inside the loop over `allResultFiles`. If no `results/single_values/??_??_50.csv` files exist, the later call to `confinement_factor_single_values(dfA, ...)` will fail because `dfA` was never defined.

Impact:
- Step 4 can crash instead of failing with a clear prerequisite error when Step 3 outputs are missing

Follow-up direction:
- Add an explicit empty-input check before concatenation and fail with a clear message

## 7. Reach-averaged geometry reconstruction in Step 5 appears broken

In `pipeline/run_confinement_values.py`, `ER_slope_margin_values(...)` calls:
- `merge_centerlines(dfSingle, _, _, False)`

Those `_` placeholders are not valid arguments in normal script execution. The call is wrapped in a broad `except`, which then falls back to `shapely.geometry.LineString()`.

Impact:
- `results/reach_averaged/*.gpkg` geometries may be empty or incorrect even when attribute calculations succeed
- The failure is silent because the exception is swallowed

Follow-up direction:
- Replace the placeholder arguments with real inputs and remove the blanket silent fallback

## 8. The DEM VRT is a required but undeclared prerequisite

Both `pipeline/main.py` and `pipeline/run_confinement_values.py` open `input_created/FAB_dem_vrt.vrt`, but I do not see active code in the current audited pipeline that writes that file. I found DEM helper code in `pipeline/dem.py` and `pipeline/select_raster.py`, but not a visible creator for this specific VRT.

Impact:
- Step 2 and Step 5 depend on filesystem state that is not yet represented as a documented active pipeline step
- A fresh run from repository code alone is not reproducible from the currently documented steps

Follow-up direction:
- Add an explicit VRT creation step or document the external prerequisite and its entrypoint clearly

## 9. Step 6 reach-averaged aggregation uses the wrong height-factor string

In `pipeline/run_confinement_values.py`, `concat_reachAveraged(directory, cross, hf)` contains:
- `if hf < 10: hf = f'02'`

That hard-codes every height factor below 10 to `02` instead of formatting the actual value.

Impact:
- Aggregated reach-averaged outputs for `hf = 0.5, 1, 1.5, 3, 4, 6, 8` will look for or write the wrong filenames
- Step 6 can silently aggregate the wrong subset of files

Follow-up direction:
- Format the actual `hf` value instead of replacing all sub-10 values with `02`

## 10. Step 6 aggregation is not robust to empty file lists

`concat_nc_conf_files(...)` and `concat_reachAveraged(...)` assume matching files exist. If none are found, they proceed to `xr.concat(dsList, ...)` or use `dfW` without guaranteeing it was assigned.

Impact:
- Step 6 can crash with unclear errors when any expected per-file Step 5 outputs are missing

Follow-up direction:
- Add explicit empty-input checks before concatenation and fail with a clear message

## 11. Step 7 only processes height factor `02`

In `pipeline/spatial_smoothing.py`, the active code hard-codes:
- `hfList = [2]`

Step 5 and Step 6 produce outputs for many height factors, but Step 7 only smooths one of them.

Impact:
- The main smoothed-results pipeline is incomplete relative to the earlier confinement-output stages
- Final clustering and final-results workflows are effectively pinned to one smoothed height factor

Follow-up direction:
- Either document `02` as the only supported smoothing height or loop over the actual intended set of height factors

## 12. Step 7 writes to `results/single_smoothed/` without visible directory creation

`pipeline/spatial_smoothing.py` writes smoothed NetCDF and pickle files to `results/single_smoothed/`, but I do not see active setup code that creates that directory.

Impact:
- Step 7 may fail on a clean filesystem

Follow-up direction:
- Add explicit directory creation before the smoothing outputs are written

## 13. Step 8 expects smoothed files for heights `02`, `03`, and `04`, but Step 7 only creates `02`

`pipeline/clustering_confinement.py` loops over `ch in [2,3,4]` and opens:
- `results/single_smoothed/global_50_0{ch}_smoothed.nc`

But the active Step 7 code only produces `global_50_02_smoothed.nc`.

Impact:
- Step 8 cannot run successfully for `03` and `04` from the currently visible upstream code path
- The clustering stage and the smoothing stage are not aligned inside the current main analysis workflow

Follow-up direction:
- Align the smoothing stage and clustering stage to the same set of height factors

## 14. Step 8 has an obvious runtime bug in the KMeans section

In `pipeline/clustering_confinement.py`, after:
- `S = KmeansTune(dfCluster, clusterCols, clusters, sampleSize, random_states)`

the code does:
- `S = [s + comb for s in S]`

but `comb` is not defined in that scope.

Impact:
- The script is likely to fail before later clustering outputs are written

Follow-up direction:
- Remove the undefined variable or replace it with the intended metadata payload
