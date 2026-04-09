#!/usr/bin/env zsh
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: ./run_world_sequential.sh <config-path> [processors] [height-factors...]"
  echo "Example: CONTINENTS='eu na sa as oc' ./run_world_sequential.sh config/paths.example_testing.json 4 2 3 4"
  exit 1
fi

CONFIG_PATH="$1"
PROCESSORS="${2:-4}"
shift 2 || true

if [[ $# -gt 0 ]]; then
  HEIGHT_FACTORS=("$@")
else
  HEIGHT_FACTORS=(2 3 4)
fi

if [[ -n "${CONTINENTS:-}" ]]; then
  CONTINENT_LIST=(${=CONTINENTS})
else
  CONTINENT_LIST=(af eu na sa as oc)
fi

echo "Running full-world canonical pipeline"
echo "  config: $CONFIG_PATH"
echo "  processors: $PROCESSORS"
echo "  continents: ${CONTINENT_LIST[*]}"
echo "  height factors: ${HEIGHT_FACTORS[*]}"

# echo
# echo "[1/8] Step 0.5 FABDEM index"
# python -m pipeline.build_fabdem_index --config "$CONFIG_PATH"

# echo
# echo "[2/8] Step 1 segmentation for all continents"
# python -m pipeline.segment_reaches \
#   "${CONTINENT_LIST[@]}" \
#   --workers 1 \
#   --config "$CONFIG_PATH"

# echo
# echo "[3/8] Step 2 canonical bend-table build for all continents"
# for CONTINENT in "${CONTINENT_LIST[@]}"; do
#   echo "  -> Step 2 continent ${CONTINENT}"
#   python -m pipeline.build_step2_results \
#     "$CONTINENT" "$PROCESSORS" \
#     --config "$CONFIG_PATH"
# done

# echo
# echo "[4/8] Step 4 confinement-factor table from all bend outputs"
# python -m pipeline.build_step4_confinement_factor \
#   --config "$CONFIG_PATH"

echo
echo "[5/8] Step 5 confinement outputs for all continents"
for HF in "${HEIGHT_FACTORS[@]}"; do
  echo "  -> height factor ${HF}"
  for CONTINENT in "${CONTINENT_LIST[@]}"; do
    echo "     -> Step 5 continent ${CONTINENT}"
    python -m pipeline.build_step5_confinement_outputs \
      "$CONTINENT" "$PROCESSORS" \
      --height-factor "$HF" \
      --config "$CONFIG_PATH"
  done
done

echo
echo "[6/8] Step 6 aggregates for all requested height factors"
for HF in "${HEIGHT_FACTORS[@]}"; do
  echo "  -> height factor ${HF}"
  python -m pipeline.build_step6_aggregates \
    --height-factor "$HF" \
    --config "$CONFIG_PATH"
done

echo
echo "[7/8] Step 7 spatial smoothing for all requested height factors"
for HF in "${HEIGHT_FACTORS[@]}"; do
  echo "  -> height factor ${HF}"
  python -m pipeline.build_step7_spatial_smoothing \
    --height-factor "$HF" \
    --workers 1 \
    --continents "${CONTINENT_LIST[@]}" \
    --config "$CONFIG_PATH"
done

echo
echo "[8/8] Step 8 clustering for all requested height factors"
python -m pipeline.build_step8_confinement_clustering \
  --height-factors "${HEIGHT_FACTORS[@]}" \
  --clusters 2 3 \
  --sample-size 10 \
  --workers 1 \
  --config "$CONFIG_PATH"

echo
echo "Optional combined Step 7-8 workflow rerun:"
echo "  python -m pipeline.build_step7_8_confinement_workflow --height-factors ${HEIGHT_FACTORS[*]} --continents ${CONTINENT_LIST[*]} --smoothing-workers 1 --clustering-workers 1 --clusters 2 3 --sample-size 10 --config $CONFIG_PATH"
echo
echo "World test finished."
