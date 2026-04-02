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

## 16. Step 2 now writes orthogonals and sampled profiles, but orchestration is still bundled in `pipeline/main.py`

Status update:
- Initial refactor implemented: `pipeline/get_orthogonals.py` now separates `build_orthogonal_lines(...)` from `sample_orthogonal_profiles(...)`
- `pipeline/main.py` now calls those as distinct substeps inside Step 2
- Step 2 now writes `results/orthogonals/*.gpkg` and `results/profiles/*.csv` before assembling `results/all/*.csv`
- A dedicated geometry-only entrypoint now exists in `pipeline/build_step2_orthogonals.py`
- The DEM-sampling stage now also exists as its own callable module in `pipeline/sample_step2_profiles.py`
- The raster-resize loop in the DEM sampling path was also fixed so it now actually increases the requested raster size when expansion is needed

Impact that remains:
- `pipeline/main.py` still owns geometry creation, DEM sampling, and final CSV assembly in one orchestration path
- The current orchestration still runs one file at a time rather than batching DEM work across saved orthogonal intermediates

Follow-up direction:
- Keep `pipeline/sample_step2_profiles.py` as the DEM stage
- Add a thin orchestration layer that can chain geometry -> profiles -> final `results/all/*.csv` without forcing them to live in one large script

## 17. Step 2 still serializes sampled profile arrays into very large CSV string columns

Status update:
- Step 3 now has a parser that can read the current nested profile encoding from `results/all/*.csv`
- That parser had to be hardened because the saved profile columns are not clean JSON; they currently include numpy-style scalar strings such as `np.float32...`

Impact that remains:
- Step 2 outputs in `results/profiles/*.csv` and `results/all/*.csv` are very large
- Step 3 remains slower and more brittle than it should be because it has to deserialize large nested list strings back into Python objects
- The current storage format makes debugging and downstream reuse harder than a structured format like Parquet or NetCDF-backed intermediates
- Simple text-line tools such as `wc -l` become misleading on downstream CSVs because nested list columns can contain embedded newlines inside quoted cells

Follow-up direction:
- Keep the current parser fix so the active pipeline still runs
- Revisit the Step 2 -> Step 3 handoff later and replace the nested-list CSV storage with a more stable intermediate format

## 18. Step 3 may be removable if Step 2 writes a bend-level structured intermediate directly

Status update:
- Step 3 currently exists mainly to transform Step 2 reach-level CSV outputs with nested list columns into one row per bend

Impact that remains:
- The pipeline carries an extra expansion stage that is expensive and tightly coupled to the current nested-list CSV format

Follow-up direction:
- Revisit the Step 2 output design later
- If Step 2 writes bend-level structured outputs directly, for example Parquet files with one row per bend, Step 3 could likely be removed or reduced to a thin compatibility layer

## 3. Step 3 now has a canonical entrypoint, but compatibility wrappers still remain

Status update:
- The canonical Step 3 entrypoint is now `pipeline/build_step3_single_values.py`
- `pipeline/open_to_single_apex.py` now delegates to that canonical CLI
- The Step 3 block in `pipeline/run_confinement_values_shell.py` now also calls the same helper instead of keeping its own file loop

Impact that remains:
- The compatibility wrapper file `pipeline/open_to_single_apex.py` still exists
- `pipeline/run_confinement_values_shell.py` still mixes Step 3 with later stages even though it no longer owns the Step 3 logic

Follow-up direction:
- Keep `pipeline/build_step3_single_values.py` as the canonical Step 3 path
- Remove or retire compatibility wrappers later when the downstream shell stages are split out

## 4. Later-stage output directories are still not created consistently outside the new Step 3 path

Status update:
- Step 3 now explicitly creates `results/single_values/` via `pipeline/build_step3_single_values.py` and `pipeline/paths.py`
- Step 5 now explicitly creates `results/reach_averaged/` and the Step 5 `results/single_values/*_conf.*` outputs via `pipeline/build_step5_confinement_outputs.py` and `pipeline/paths.py`

Impact that remains:
- Later stages still write to directories such as `results/single_smoothed/` without one consistent setup path across the active pipeline

Follow-up direction:
- Continue the same path/setup cleanup pattern for Step 6 and later stages

## 5. Step 4 now has a canonical entrypoint, but later stages still consume it through older directory-based code

Status update:
- The canonical Step 4 entrypoint is now `pipeline/build_step4_confinement_factor.py`
- The Step 4 block in `pipeline/run_confinement_values_shell.py` now calls that helper instead of keeping its own concatenation logic

Impact that remains:
- The shell script still mixes Step 4 with later stages even though it no longer owns the Step 4 logic
- The broader final workflow is still spread across multiple stage entrypoints plus the older shell wrapper

Follow-up direction:
- Keep `pipeline/build_step4_confinement_factor.py` as the canonical Step 4 path
- Continue splitting Step 5 and later stages into their own explicit entrypoints and orchestration layers

## 6. Step 4 is not robust to an empty Step 3 result set

Status update:
- The canonical Step 4 helper now checks for missing inputs explicitly
- When no Step 3 files are available, it raises a clear `FileNotFoundError` instead of failing later with an undefined dataframe

Impact that remains:
- No specific empty-input bug remains in the canonical Step 4 path
- Broader downstream robustness work is still needed in later stages

Follow-up direction:
- Keep the explicit prerequisite checks as later stages are split out

## 7. Reach-averaged geometry reconstruction in Step 5 appears broken

In `pipeline/run_confinement_values.py`, `ER_slope_margin_values(...)` calls:
- `merge_centerlines(dfSingle, _, _, False)`

Those `_` placeholders are not valid arguments in normal script execution. The call is wrapped in a broad `except`, which then falls back to `shapely.geometry.LineString()`.

Impact:
- `results/reach_averaged/*.gpkg` geometries may be empty or incorrect even when attribute calculations succeed
- The failure is silent because the exception is swallowed

Follow-up direction:
- Replace the placeholder arguments with real inputs and remove the blanket silent fallback

## 8. FABDEM source indexing is now explicit, but later DEM-dependent code still needs broader path cleanup

Status update:
- Initial fix implemented: `pipeline/build_fabdem_index.py` now creates `input_created/dem/FAB_dem_vrt.vrt` and `input_created/dem/FAB_dem_bounds.gpkg`
- The source FABDEM tile directory is now intended to be configured once via `fabdem_dir`

Remaining follow-up:
- The broader DEM-related workflow still mixes centralized FABDEM paths with older hard-coded base-directory patterns in other active scripts

Impact that remains:
- The VRT is no longer undeclared, but the rest of the pipeline has not yet been fully migrated to one consistent path/config model

## 9. Step 6 height-factor formatting is fixed in the canonical aggregation path

Status update:
- `concat_reachAveraged(...)` now uses the same height-factor formatting helper as Step 5
- The canonical Step 6 entrypoint is now `pipeline/build_step6_aggregates.py`

Impact that remains:
- No specific height-factor formatting bug remains in the canonical Step 6 path
- The older shell wrapper still exists, but it now calls the canonical Step 6 helper

Follow-up direction:
- Keep the shared height-factor formatting logic aligned across Steps 5 and 6

## 10. Step 6 empty-input handling is fixed in the canonical aggregation path

Status update:
- `concat_nc_conf_files(...)` and `concat_reachAveraged(...)` now check for missing Step 5 inputs explicitly
- The canonical Step 6 path now raises clear `FileNotFoundError` messages when no matching per-file outputs exist for the requested height factor

Impact that remains:
- No specific empty-input bug remains in the canonical Step 6 path
- Broader later-stage robustness work is still needed in Step 7 and beyond

Follow-up direction:
- Keep the explicit prerequisite checks as the remaining downstream stages are split out

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
