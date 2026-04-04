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

# end of issues
