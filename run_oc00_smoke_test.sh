#!/usr/bin/env zsh
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: ./run_oc00_smoke_test.sh <config-path>"
  echo "Example: ./run_oc00_smoke_test.sh config/paths.local.oc00_test.json"
  exit 1
fi

CONFIG_PATH="$1"
VECTOR_FILE="results/new_segments/vector/oc_00_reach_new_segments.gpkg"
BEND_FILE="results/bends/oc_00_50.parquet"

echo
echo "[1/8] Step 1.5 FABDEM index"
python -m pipeline.build_fabdem_index --config "$CONFIG_PATH"

echo
echo "[2/8] Step 1 segmentation for oc"
python -m pipeline.segment_reaches oc --workers 1 --config "$CONFIG_PATH"

echo
echo "[3/8] Step 2 canonical bend-table build for oc_00"
python -m pipeline.build_step2_results \
  --vector-file "$VECTOR_FILE" \
  --config "$CONFIG_PATH"

echo
echo "[4/8] Step 4 confinement-factor table from bend-level Step 2 output"
python -m pipeline.build_step4_confinement_factor \
  --input-file "$BEND_FILE" \
  --config "$CONFIG_PATH"

echo
echo "[5/8] Step 5 confinement outputs for height factor 2"
python -m pipeline.build_step5_confinement_outputs \
  --input-file "$BEND_FILE" \
  --height-factor 2 \
  --config "$CONFIG_PATH"

echo
echo "[6/8] Step 6 aggregates for height factor 2"
python -m pipeline.build_step6_aggregates \
  --height-factor 2 \
  --config "$CONFIG_PATH"

echo
echo "[7/8] Step 7 spatial smoothing for height factor 2"
python -m pipeline.build_step7_spatial_smoothing \
  --height-factor 2 \
  --workers 1 \
  --continents oc \
  --config "$CONFIG_PATH"

echo
echo "[8/8] Step 8 clustering smoke test for height factor 2"
python -m pipeline.build_step8_confinement_clustering \
  --height-factors 2 \
  --clusters 2 3 \
  --sample-size 10 \
  --workers 1 \
  --config "$CONFIG_PATH"

echo
echo "Smoke test finished."
echo "Optional compatibility export:"
echo "  python -m pipeline.build_step3_single_values --input-file $BEND_FILE --config $CONFIG_PATH"
echo
echo "Optional full Step 7-8 workflow for 2 3 4 after building Step 5 and 6 for all three height factors:"
echo "  python -m pipeline.build_step7_8_confinement_workflow --height-factors 2 3 4 --continents oc --smoothing-workers 1 --clustering-workers 1 --clusters 2 3 --sample-size 10 --config $CONFIG_PATH"
