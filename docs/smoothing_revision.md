# Step 7 smoothing revision

The existing `bend_neighbor_graph` and `smooth_attributes` functions are retained
unchanged for baseline reproduction. Step 7 defaults to `--method legacy` until
sensitivity results justify adopting a new setting. No final N or alpha has been
chosen. Steps 2–6 are not modified.

## New method

`pipeline/local_bend_smoothing.py` orders bends by descending `bendDistOut` within
continent/network/reach. Rank 1 is upstream. Interior links follow adjacent ranks;
reach boundaries follow the existing `combined_reach_up`/`combined_reach_dn`
references. At confluences this follows the upstream reach already selected by the
pipeline, not every incoming tributary. Missing or out-of-network references stop
the walk and are counted. Repeated nodes fail explicitly for inspection. Invalid
lengths, inconsistent reach references, and ambiguous distance ordering also fail.

The averaging set contains the focal bend once, plus at most N upstream and N
downstream bends. Missing neighbors on one side are not replaced on the other.
No 20 km selection cutoff applies to this method.

For focal length L, each bend on a path uses effective length `max(bendLen, L)`
when the length floor is enabled. Adjacent effective half-lengths are summed along
the path; these are weighting distances, not actual geographic midpoint distances.
The bandwidth is `alpha * L`. Weights are `exp(-0.5 * (distance / bandwidth)**2)`
and normalized to sum to one, including the focal bend at distance zero.
The existing weighted mean, NaN propagation, singleton zero STD, and sample-count
STD correction are preserved. Infinite attributes fail rather than silently enter
clustering. STD remains sensitive to candidate count even for tiny weights.

This floor limits short-neighbor influence; it does not make pairwise influence
reciprocal. Long neighbors keep their true length. The three tunable choices are
candidate count N, bandwidth alpha, and whether the floor is enabled.

## Commands

Run with an environment containing the repository's scientific dependencies. On
the current workstation `/opt/anaconda3/envs/test-orthogonals/bin/python` was used.
Set `PYTHONDONTWRITEBYTECODE=1` to avoid modifying tracked bytecode files.

Single experimental Step 7 run (settings are provisional):

```sh
python -m pipeline.build_step7_spatial_smoothing \
  --config /path/to/paths.json --height-factor 2 --continents oc \
  --method local --neighbors-per-direction 3 --alpha 0.75 --length-floor \
  --output-dir results/smoothing_revision/example/single_smoothed
```

The combined Step 7–8 wrapper accepts the same method parameters and an experiment
root `--output-dir`; it places smoothing files under `single_smoothed` and tuning
scores at the experiment root. Step 8 also accepts `--input-dir` and `--output-dir`.
Existing tuning behavior is preserved. Step 7 concatenates only continents
produced by that invocation. A settings manifest rejects reuse of a directory for
a different smoothing configuration.

Nine-setting grid, with an existing baseline:

```sh
python -m pipeline.build_smoothing_sensitivity \
  --input /path/to/global_50_02_conf.nc \
  --baseline /path/to/verified_baseline.nc \
  --output-dir results/smoothing_revision/grid_hf2 \
  --neighbors 2 3 4 --alphas 0.5 0.75 1
```

The output directory must be new. Omit `--continents` for all input continents.
Each invocation handles one supplied height-factor input; repeat for other required
height factors. Add `--run-legacy` instead of `--baseline` to recompute the exact
current method and write an audit of selected real bend IDs and weights. This can
be expensive on large networks. A pre-existing archive is not proof of how it was
generated; check provenance and matched IDs before calling it the production baseline.

Add `--clusters <k values>` to retain assignments and compare clustering. Fixed
settings: K-means n_init=10; GMM tied covariance, n_init=5, max_iter=200, kmeans
initialization. Both use the requested seeds (default 20,43,50). These are explicit
comparison settings, not a claim to reproduce the manuscript's final model.
Use the relevant k values once identified. Existing eight-attribute scaling and
three-component PCA preparation are refit for each configuration on an identical,
finite common bend population. Score samples use identical bend IDs for each seed.
Eligibility counts are reported before intersection. Labels are aligned by maximum
overlap for proportions; ARI needs no alignment. Results include non-convergence.

After selecting a candidate, add `--no-floor-control N ALPHA` to isolate the floor's
effect. A pilot may include a provisional control, but this does not select it as
the final method. `--max-network-bends` and `--network-limit` are diagnostic subset
options; such results must not be reported as continent/global robustness.

Outputs include a manifest, network sizes, topology diagnostics, per-configuration
NetCDF, slope summaries, all pairwise slope comparisons (global, continent, input
length quartiles), separate shared/per-setting timings, and optional assignment,
score, ARI, eligibility, and aligned cluster-proportion CSV files.

Full-continent benchmark in independent fresh processes:

```sh
python -m pipeline.benchmark_bend_smoothing \
  --input /path/to/global_50_02_conf.nc --continents oc \
  --output-dir results/smoothing_revision/benchmark \
  --neighbors 3 --alpha 0.75 --legacy-timeout 120 --legacy-memory-gb 4
```

Reports shared loading/preparation separately from graph construction, averaging,
and output writing. Peak process RSS includes the loaded input and Python runtime.
Legacy timeout/memory-limit results are explicitly incomplete; they are not actual
completed before/after speed ratios. Increase limits only if resources permit.

## Verification

```sh
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -v
```

Tests cover analytic equal-length and short-neighbor weights, longer neighbors,
network/reach boundaries, original arithmetic parity, singleton/NaN behavior,
invalid lengths, and repeated-node detection. The pilot and benchmark artifacts in
`results/smoothing_revision` document actual runs separately from synthetic tests.
