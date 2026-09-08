# SWORD–VBET validation pilot

This workflow compares the manuscript's SWORD/FABDEM confinement result with
independently produced VBET valley-bottom widths. It deliberately begins from a
Riverscapes **RME IGO raw export**, not VBET valley-bottom polygons.

## Scope

Use the 2025 CONUS run only. It covers the lower 48 where model runs completed;
it does not provide systematic coverage for Canada, Mexico, Alaska, or
transboundary HUCs. Coverage is established by the returned IGO points, not by
assuming every SWORD reach in a bounding box is covered.

## One-time pilot export

In Riverscapes Reports, draw a small pilot AOI (one HUC2 or several adjacent
HUC8s) and request an IGO/raw data export. Before a large export, confirm it
contains:

- a unique IGO identifier and point geometry;
- full-valley-bottom `integrated_width` (metres), or its documented current
  RME field name;
- optional low-lying floodplain width and centreline-length fields;
- any stream-size/drainage-area fields useful for later filtering.

Save the export as a GeoPackage. Do not download VBET DGO polygons for the
primary analysis; they are only needed for visual quality assurance.

## Commands

Use the names present in your own SWORD and RME export. The illustrative names
below must be replaced after inspecting their schemas.

```bash
# Filter only valid manuscript records; this does not claim VBET coverage.
python -m pipeline.build_vbet_validation candidates \
  --sword results/reach_averaged/na_reaches.gpkg \
  --output results/vbet_validation/sword_candidates.gpkg \
  --reach-id reach_id \
  --metric confinement_ratio

# Match IGO points by geometry. Start at 250 m, then inspect match distances.
python -m pipeline.build_vbet_validation match \
  --sword results/vbet_validation/sword_candidates.gpkg \
  --igos data/vbet/pilot_igos.gpkg \
  --output results/vbet_validation/pilot_igo_matches.gpkg \
  --reach-id reach_id \
  --igo-id igo_id \
  --max-distance-m 250

# Construct reach-level VBET/SWORD width ratios.
python -m pipeline.build_vbet_validation aggregate \
  --matches results/vbet_validation/pilot_igo_matches.gpkg \
  --output results/vbet_validation/pilot_reach_widths.parquet \
  --reach-id reach_id \
  --width integrated_width \
  --denominator bend_width_m \
  --min-igos 3
```

The output includes median VBET width, its IQR, the median IGO-to-SWORD match
distance, the number of IGOs, and the VBET/SWORD width ratio. It uses the median
IGO width per reach for the pilot. This is robust to isolated broad valley
segments; a length-weighted aggregation can be added after we confirm the
specific RME metric fields and whether the paper's final denominator is
reach-level or bend-level.

## Quality gates before expansion

1. Visually inspect 100--200 matches, balanced across confinement classes,
   stream sizes, and regions.
2. Plot and set a justified threshold from the observed match-distance
   distribution; do not retain ambiguous confluence/distributary matches.
3. Check that at least three IGOs support each retained reach and that no small
   set of large rivers dominates the sample.
4. Compare the continuous width ratios and the manuscript's confinement classes
   only after passing these spatial checks.

For the scale-up, export IGOs by HUC2/HUC4 batches and retain the source-export
metadata, query settings, and match-summary CSV alongside each batch.
