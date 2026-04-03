# Known Issues

## 1. Hard-coded base paths

Status update:
- The canonical Steps 1 through 8 now run through `pipeline/paths.py` and `config/paths.local.json`
- `pipeline/run_confinement_values_shell.py`, `pipeline/select_raster.py`, and the active path-setup cells in `pipeline/final_results.ipynb` no longer hard-code `/scratch/6256481/` or `/Volumes/...`

Impact that remains:
- The main pipeline is now portable through the shared path/config layer
- Some older helper functions still carry compatibility-style `directory` arguments or legacy path assumptions outside the canonical stage flow
- A few notebook comments still mention old absolute example paths

Examples that still need later cleanup:
- `pipeline/dem.py` MERIT-specific helper paths
- older compatibility-oriented helpers inside `pipeline/run_confinement_values.py`
- stale commented path examples in `pipeline/final_results.ipynb`

## 2. Step 1 entrypoint separation

Status update:
- Initial fix implemented: Step 1 now has a dedicated entrypoint in `pipeline/segment_reaches.py`
- Step 1 paths are now loaded through `pipeline/paths.py` and `config/paths.local.json`

Remaining follow-up:
- Keep the canonical stage entrypoints on the shared path system and continue pruning older helper/compatibility code around them

Impact that remains:
- Step 1 is now cleanly separated and the canonical later stages now use the same project layout
- Remaining work is mostly around older helper modules and convenience wrappers, not the main audited pipeline path

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
- Step 7 now explicitly creates `results/single_smoothed/` via `pipeline/spatial_smoothing.py` and `pipeline/paths.py`
- Step 8 now writes its score tables through the shared results root in `pipeline/clustering_confinement.py`

Impact that remains:
- No specific later-stage output-directory bug remains in the canonical Steps 3 through 8 path
- The same cleanup pattern is still needed in older non-canonical helpers and notebooks

Follow-up direction:
- Keep using `pipeline/paths.py` as older helper scripts are either cleaned or retired

## 5. Step 4 now has a canonical entrypoint, but the multi-stage shell wrapper still exists as compatibility code

Status update:
- The canonical Step 4 entrypoint is now `pipeline/build_step4_confinement_factor.py`
- The Step 4 block in `pipeline/run_confinement_values_shell.py` now calls that helper instead of keeping its own concatenation logic
- `pipeline/run_confinement_values_shell.py` is now path-aware and delegates to the canonical Step 3 to Step 6 entrypoints

Impact that remains:
- The shell script still mixes several later stages in one compatibility wrapper
- The broader final workflow is still spread across multiple explicit stage entrypoints plus that wrapper

Follow-up direction:
- Keep `pipeline/build_step4_confinement_factor.py` as the canonical Step 4 path
- Retire or minimize the compatibility wrapper later once a preferred high-level orchestration path is settled

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

## 11. Step 7 height-factor handling is fixed in the canonical smoothing path

Status update:
- The canonical Step 7 entrypoint is now `pipeline/build_step7_spatial_smoothing.py`
- The smoothing stage now accepts `--height-factor` instead of hard-coding `hfList = [2]`

Impact that remains:
- No specific hard-coded `02` bug remains in the canonical Step 7 path
- There is still no higher-level workflow wrapper that automatically runs the full intended set of smoothed height factors for downstream analysis

Follow-up direction:
- Decide which height factors the final workflow actually requires, then orchestrate those explicitly

## 12. Step 7 output-directory creation is fixed in the canonical smoothing path

Status update:
- The canonical Step 7 path now creates `results/single_smoothed/` before writing smoothed NetCDF and pickle outputs

Impact that remains:
- No specific clean-filesystem directory bug remains in the canonical Step 7 path
- The older downstream scripts still need the same path/setup review pattern

Follow-up direction:
- Keep using `pipeline/paths.py` for later-stage output setup instead of ad hoc directory strings

## 13. Step 8 expects smoothed files for heights `02`, `03`, and `04`, but there is no matching multi-height Step 7 orchestration yet

Status update:
- The canonical Step 8 entrypoint is now `pipeline/build_step8_confinement_clustering.py`
- The clustering stage now accepts `--height-factors` instead of hard-coding `ch in [2,3,4]`
- The canonical Step 8 path reads `results/single_smoothed/global_50_{hf}_smoothed.nc` through the shared path setup

The canonical Step 7 path can now smooth any requested height factor, and the canonical Step 8 path can now cluster any requested height factor, but I do not see orchestration in the current repo that automatically produces `02`, `03`, and `04` together before clustering runs.

Impact:
- Step 8 cannot run successfully for `03` and `04` from the currently visible upstream code path
- The clustering stage and the smoothing stage are not aligned inside the current main analysis workflow

Follow-up direction:
- Align the smoothing stage and clustering stage to the same set of height factors

## 14. Step 8 KMeans runtime bug is fixed in the canonical clustering path

Status update:
- The canonical Step 8 path now writes KMeans score tables directly from `KmeansTune(...)`
- The undefined `comb` variable has been removed from the canonical clustering code

Impact that remains:
- No specific `comb` runtime bug remains in the canonical Step 8 path
- Broader clustering-workflow cleanup is still needed where older notebooks and scripts expect the previous hard-coded paths

Follow-up direction:
- Keep the canonical Step 8 path as the source of truth and update any downstream notebooks to its current paths and outputs


# end of issues
