# Resolved Issues

These issues are considered resolved in the canonical pipeline path. Some older compatibility helpers or notebook comments may still contain historical patterns, but the audited Step 1 to Step 8 workflow no longer depends on the original broken behavior.

## 1. Hard-coded base paths

Status update:
- The canonical Steps 1 through 8 now run through `pipeline/paths.py` and `config/paths.local.json`
- the former Step 3 to Step 6 shell wrapper, `pipeline/select_raster.py`, and the active path-setup cells in `pipeline/final_results.ipynb` no longer hard-code `/scratch/6256481/` or `/Volumes/...`

Resolution:
- The main confinement pipeline is now portable through the shared path/config layer
- Remaining MERIT-specific helper paths and compatibility-style `directory` arguments do not affect the current active workflow

## 2. Step 1 entrypoint separation

Status update:
- Step 1 now has a dedicated entrypoint in `pipeline/segment_reaches.py`
- Step 1 paths are now loaded through `pipeline/paths.py` and `config/paths.local.json`
- `pipeline/main.py` no longer owns segmentation and now only delegates to the canonical Step 2 entrypoint

Resolution:
- The Step 1 entrypoint separation problem is fixed in the canonical pipeline path

## 16. Step 2 orchestration split out of `pipeline/main.py`

Status update:
- The canonical Step 2 entrypoint is now `pipeline/build_step2_results.py`
- The active Step 2 worker/orchestration logic now lives in `pipeline/step2.py`
- `pipeline/build_step2_orthogonals.py` remains the geometry-only entrypoint
- `pipeline/sample_step2_profiles.py` remains the standalone DEM/profile sampling stage
- `pipeline/main.py` is now only a compatibility wrapper that delegates to the canonical Step 2 entrypoint

Resolution:
- The canonical Step 2 path no longer bundles orchestration in `pipeline/main.py`

## 3. Step 3 compatibility wrappers

Status update:
- The canonical Step 3 entrypoint is `pipeline/build_step3_single_values.py`
- The old Step 3 wrapper files are no longer part of the active `pipeline/` tree

Resolution:
- The canonical Step 3 path no longer depends on compatibility wrappers inside the active pipeline module

## 4. Later-stage output-directory creation

Status update:
- Step 3 now explicitly creates `results/single_values/` via `pipeline/build_step3_single_values.py` and `pipeline/paths.py`
- Step 5 now explicitly creates `results/reach_averaged/` and the Step 5 `results/single_values/*_conf.*` outputs via `pipeline/build_step5_confinement_outputs.py` and `pipeline/paths.py`
- Step 7 now explicitly creates `results/single_smoothed/` via `pipeline/spatial_smoothing.py` and `pipeline/paths.py`
- Step 8 now writes its score tables through the shared results root in `pipeline/clustering_confinement.py`

Resolution:
- No specific later-stage output-directory bug remains in the canonical Steps 3 through 8 path

## 6. Step 4 empty-input handling

Status update:
- The canonical Step 4 helper now checks for missing inputs explicitly
- When no canonical bend-level Step 2 files are available, it raises a clear `FileNotFoundError` instead of failing later with an undefined dataframe

Resolution:
- No specific empty-input bug remains in the canonical Step 4 path

## 5. Step 4 compatibility wrapper duplication

Status update:
- The canonical Step 4 entrypoint is `pipeline/build_step4_confinement_factor.py`
- The old multi-stage shell wrapper has been moved out of the active `pipeline/` tree

Resolution:
- The canonical Step 4 path no longer depends on the old compatibility wrapper inside the active pipeline module

## 8. FABDEM source indexing and active DEM path cleanup

Status update:
- `pipeline/build_fabdem_index.py` creates `input_created/dem/FAB_dem_vrt.vrt` and `input_created/dem/FAB_dem_bounds.gpkg`
- Active FABDEM-based code now reads those shared derived inputs through `pipeline/paths.py`
- The old MERIT-specific helper functions were moved out of `pipeline/dem.py` into `old/merit_dem.py`

Resolution:
- No active DEM-path ambiguity remains in the canonical FABDEM-based workflow

## 9. Step 6 height-factor formatting

Status update:
- `concat_reachAveraged(...)` now uses the same height-factor formatting helper as Step 5
- The canonical Step 6 entrypoint is now `pipeline/build_step6_aggregates.py`

Resolution:
- No specific height-factor formatting bug remains in the canonical Step 6 path

## 10. Step 6 empty-input handling

Status update:
- `concat_nc_conf_files(...)` and `concat_reachAveraged(...)` now check for missing Step 5 inputs explicitly
- The canonical Step 6 path now raises clear `FileNotFoundError` messages when no matching per-file outputs exist for the requested height factor

Resolution:
- No specific empty-input bug remains in the canonical Step 6 path

## 11. Step 7 height-factor handling

Status update:
- The canonical Step 7 entrypoint is now `pipeline/build_step7_spatial_smoothing.py`
- The smoothing stage now accepts `--height-factor` instead of hard-coding `hfList = [2]`

Resolution:
- No specific hard-coded `02` bug remains in the canonical Step 7 path

## 12. Step 7 output-directory creation

Status update:
- The canonical Step 7 path now creates `results/single_smoothed/` before writing smoothed NetCDF and pickle outputs

Resolution:
- No specific clean-filesystem directory bug remains in the canonical Step 7 path

## 14. Step 8 KMeans runtime bug

Status update:
- The canonical Step 8 path now writes KMeans score tables directly from `KmeansTune(...)`
- The undefined `comb` variable has been removed from the canonical clustering code

Resolution:
- No specific `comb` runtime bug remains in the canonical Step 8 path

## 19. Redundant compatibility interfaces and helper cleanup

Status update:
- The old wrapper files `pipeline/open_to_single_apex.py` and `pipeline/run_confinement_values_shell.py` were moved into `old/`
- The old MERIT-only helpers were moved from `pipeline/dem.py` into `old/merit_dem.py`
- The active Step 5 and Step 6 code in `pipeline/run_confinement_values.py` now only uses the canonical config-driven path layer

Resolution:
- The main active pipeline no longer carries the redundant compatibility interfaces that were tracked in this issue

## 15. Step 1 multi-file-per-continent guard wording

Status update:
- The canonical Step 1 path intentionally validates that there is exactly one reach file and one node file per continent
- The validation error now states clearly that the configured input folders contain too many files for that continent and lists the matching filenames

Resolution:
- The one-file-per-continent check is now treated as an explicit guardrail, not as an ambiguous pipeline failure

## 7. Reach-averaged geometry reconstruction in Step 5

Status update:
- The broken placeholder `merge_centerlines(dfSingle, _, _, False)` call has been removed from the canonical Step 5 path
- The canonical bend-level Step 2 output now carries `centerlineWkt`
- Step 5 now uses that saved centerline to populate reach-averaged geometry in the canonical path

Resolution:
- The canonical Step 5 path no longer depends on the broken placeholder merge call for reach-averaged geometry

## 17. Step 2 nested profile CSV serialization

Status update:
- The canonical Step 2 path now writes `results/bends/*.parquet` with one row per bend
- The DEM profile arrays now travel through the canonical pipeline as bend-level columns instead of large reach-level CSV string blobs
- `pipeline/sample_step2_profiles.py` still supports CSV export for debugging, but that is no longer the canonical handoff

Resolution:
- The canonical pipeline no longer depends on the large nested-list `results/profiles/*.csv` and `results/all/*.csv` handoff format

## 18. Step 3 removable after bend-level Step 2 output

Status update:
- The canonical Step 2 path now writes bend-level tables directly
- Step 4 and Step 5 now consume those bend-level Step 2 outputs directly
- `pipeline/build_step3_single_values.py` now exists only as a compatibility exporter to the legacy `results/single_values/*.csv` format
- The canonical bend-level builder now creates one row per bend directly instead of flattening reach rows with `explode(...)`

Resolution:
- The canonical pipeline no longer requires Step 3 as a mandatory transformation stage

## 13. Step 7 and Step 8 multi-height orchestration

Status update:
- The canonical Step 7 entrypoint can smooth any requested height factor
- The canonical Step 8 entrypoint can cluster any requested height factor
- The canonical workflow `pipeline/build_step7_8_confinement_workflow.py` now runs the required height-factor set together, by default `2 3 4`

Resolution:
- The canonical analysis workflow no longer relies on manual alignment between Step 7 smoothing and Step 8 clustering
