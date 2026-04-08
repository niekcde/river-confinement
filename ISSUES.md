# Known Issues

Resolved issues now live in [RESOLVED_ISSUES.md](/Users/6256481/Code/river-confinement/RESOLVED_ISSUES.md).

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
- The current pipeline still writes many large intermediate artifacts across `results/new_segments/`, `results/orthogonals/`, `results/bends/`, optional compatibility exports in `results/single_values/`, and later outputs
- Some of those intermediates still exist mainly to support debugging or compatibility rather than the now-simpler canonical handoff

Expected impact:
- Lower disk usage
- Less duplicated IO
- Simpler reruns once the preferred pipeline boundaries are settled

Follow-up direction:
- Decide which intermediates are truly required for reproducibility/debugging
- Remove or collapse the rest after the remaining functional issues are addressed
- Revisit this together with the remaining Step 5 and final-workflow issues

## 22. Step 2 geometry stage is still slow

Status update:
- The canonical Step 2 path now works and writes bend-level Parquet, but the geometry build remains expensive
- In the `oc_00` smoke test, the canonical Step 2 run still took roughly half an hour for one file

Why this matters:
- Step 2 remains the dominant runtime bottleneck before any full-continent or global rerun
- The current geometry stage still does substantial per-reach work in Python loops:
  - centerline merging
  - smoothing
  - inflection-point extraction
  - orthogonal generation

Follow-up direction:
- Profile the Step 2 geometry path directly
- Identify the hottest functions in `pipeline/step2.py`, `pipeline/inflection_points.py`, `pipeline/smoothing.py`, and `pipeline/get_orthogonals.py`
- Reduce repeated geometry conversions and per-row dataframe work where possible

## 23. Step 7 smoothing is still computationally heavy and loop-based

Status update:
- The canonical Step 7 path works, but smoothing still builds neighbor graphs and applies smoothing in explicit Python loops
- Even the small `oc_00` smoke test spent noticeable time in graph construction and per-bend smoothing

Why this matters:
- Full-continent or global reruns will scale this cost up substantially
- The current implementation in `pipeline/spatial_smoothing.py` depends on repeated dataframe filtering and iterative assignment

Follow-up direction:
- Profile `bend_neighbor_graph(...)` and `smooth_attributes(...)`
- Reduce repeated dataframe lookups keyed by `bendID`
- Consider a more compact graph/smoothing representation once the final feature set is stable

## 24. Step 2 smoothing fallback warnings need review

Status update:
- The smoke test emitted repeated messages like `Smoothing broken at window size 2 (...) (...)`
- The run still completed successfully, so this is not currently a blocker

Why this matters:
- It is still unclear whether those warnings indicate expected edge cases or a real geometry-quality problem
- If they represent recoverable failures, the pipeline should classify or summarize them more clearly

Follow-up direction:
- Inspect the warning path in `pipeline/smoothing.py`
- Decide whether these cases are acceptable fallback behavior or should be treated as a data-quality/debug issue

## 25. Step 6 still emits an xarray concat future-warning

Status update:
- The canonical Step 6 run works, but the smoke test emitted an `xarray` `FutureWarning` during `xr.concat(...)`

Why this matters:
- The current behavior is stable now, but future `xarray` releases may change default concat behavior
- This is low risk today, but easy to make explicit now

Follow-up direction:
- Update `concat_nc_conf_files(...)` to pass explicit concat options instead of relying on defaults

## 26. The canonical pipeline has only been smoke-tested on `oc_00`

Status update:
- The smallest-file smoke test now runs through Steps 1.5 to 8 successfully
- The audited output structure looks correct for `oc_00`

Why this matters:
- This is still not equivalent to a continent-scale or global validation
- Some batching, memory, or aggregation issues may only appear when many files are processed together

Follow-up direction:
- Run a larger validation pass, for example all `oc` files and then one broader multi-continent run
- Confirm that the canonical workflow also behaves correctly for height factors `2 3 4`

## 27. Step 8 writes score tables, but not a canonical cluster-assignment output

Status update:
- The canonical Step 8 path writes:
  - `KmeanScores_confinement_{hf}.csv`
  - `GMMSCORES_confinement_{hf}.csv`
- The notebook still performs the actual clustering assignments interactively from the smoothed dataset

Why this matters:
- If the final workflow should be reproducible outside the notebook, cluster labels and probabilities likely need a canonical saved output
- Right now Step 8 tunes/evaluates clustering, but does not materialize a standard downstream cluster dataset

Follow-up direction:
- Decide whether the canonical pipeline should export:
  - one clustered NetCDF/Parquet table
  - or one per-height-factor cluster-assignment file
- Align that decision with how `pipeline/final_results.ipynb` uses clustering outputs

## 28. `final_results.ipynb` still contains a lot of manual analysis logic

Status update:
- The notebook now points at the canonical path setup and smoothed outputs
- The main stale input assumptions were cleaned up
- But the notebook still contains a large amount of exploratory/manual code mixed with final figure generation

Why this matters:
- The final workflow is not yet fully explicit or minimal
- Some figure-generation steps still depend on manual execution order and interactive clustering steps

Follow-up direction:
- Decide which notebook cells are canonical final-analysis steps and which are exploratory leftovers
- Reduce or separate exploratory analysis from the final figure workflow
- If needed, promote the truly canonical final-analysis steps into dedicated scripts

## 29. The smoke-test helper script is linear and does not support resumable stages

Status update:
- `run_oc00_smoke_test.sh` now works for the canonical smallest-file test
- It currently runs the stages in a simple linear sequence

Why this matters:
- Debugging often requires restarting from a later successful stage
- The current script is good for a clean smoke test, but not yet ideal for repeated development testing

Follow-up direction:
- Add optional stage skipping or start-from-stage flags if repeated smoke testing becomes common
- Keep the script small and focused on the canonical path

## 30. Step 4 lookup-table size could probably be reduced

Status update:
- The canonical Step 4 path works and writes `confinement_factor_50.csv`
- The current table still contains one row per bend-width instance from the bend dataset

Why this matters:
- The lookup table is conceptually a width-to-factor mapping, not necessarily a full bend-level export
- A deduplicated or binned representation may reduce file size and simplify reuse

Follow-up direction:
- Decide whether Step 4 should continue storing the full per-bend width series
- Or switch to a deduplicated lookup representation if it does not change Step 5 behavior

# end of issues
