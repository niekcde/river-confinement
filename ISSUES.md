# Known Issues

Resolved issues now live in [RESOLVED_ISSUES.md](/Users/6256481/Code/river-confinement/RESOLVED_ISSUES.md).

## 15. Step 1 still assumes one reach file and one node file per continent

The new `pipeline/segment_reaches.py` entrypoint currently fails if more than one `*17*.gpkg` reach file or node file is found for a continent.

Why this exists:
- The downstream active pipeline still assumes Step 1 outputs are named like `{continent}_{group}_...`
- That naming scheme comes from `new_reach_definition(...)`, where `group` is the segmented network-group id used later by `pipeline/main.py`

Impact:
- Step 1 is now stable and explicit, but it is not yet generalized to multi-file-per-continent source layouts
- Supporting multiple source files per continent will require aligning the Step 2 file-batching logic and the later file-sorting logic, not just Step 1

Follow-up direction:
- Keep the current one-file-per-continent assumption for now
- Revisit this only when the later stages are refactored onto the same path/file naming system

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

## 7. Reach-averaged geometry reconstruction in Step 5 appears broken

In `pipeline/run_confinement_values.py`, `ER_slope_margin_values(...)` calls:
- `merge_centerlines(dfSingle, _, _, False)`

Those `_` placeholders are not valid arguments in normal script execution. The call is wrapped in a broad `except`, which then falls back to `shapely.geometry.LineString()`.

Impact:
- `results/reach_averaged/*.gpkg` geometries may be empty or incorrect even when attribute calculations succeed
- The failure is silent because the exception is swallowed

Follow-up direction:
- Replace the placeholder arguments with real inputs and remove the blanket silent fallback

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

## 20. Replace segmented Step 1 source discovery with one global Parquet input

Status update:
- This is a planned future input-model change, not a current pipeline bug

Planned direction:
- Move away from continent-by-continent `*17*.gpkg` source discovery as the main long-term input mode
- Introduce one large structured Parquet input file as the canonical upstream source
- Treat that Parquet file as the single Step 1 input instead of discovering multiple vector/node source files at runtime

Expected impact:
- Simpler Step 1 orchestration
- Less fragile file discovery logic
- Easier later batching and partitioning strategy

Follow-up direction:
- Revisit Step 1 and Step 2 together when the new unified Parquet source is ready

## 21. Reduce the number of intermediate saved files across the pipeline

Status update:
- This is a planned cleanup/performance direction, not a single current runtime bug

Why this matters:
- The current pipeline writes many large intermediate artifacts across `results/new_segments/`, `results/orthogonals/`, `results/profiles/`, `results/all/`, `results/single_values/`, and later outputs
- Some of those intermediates exist mainly to support the current staged refactor or current serialization format

Expected impact:
- Lower disk usage
- Less duplicated IO
- Simpler reruns once the preferred pipeline boundaries are settled

Follow-up direction:
- Decide which intermediates are truly required for reproducibility/debugging
- Remove or collapse the rest after the remaining functional issues are addressed
- Revisit this together with Issue 17 and Issue 18

# end of issues
