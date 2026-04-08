#!/usr/bin/env zsh
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: ./run_continent_test.sh <config-path> [continent] [processors] [height-factor]"
  echo "Example: ./run_continent_test.sh config/paths.local.af_test.json af 4 2"
  exit 1
fi

CONFIG_PATH="$1"
CONTINENT="${2:-af}"
PROCESSORS="${3:-1}"
HEIGHT_FACTOR="${4:-2}"

python - "$CONFIG_PATH" "$CONTINENT" <<'PY'
import sys

from pipeline.paths import load_project_paths

config_path = sys.argv[1]
target_continent = sys.argv[2]
paths = load_project_paths(config_path)

print(f"Using results_root: {paths.results_root}")

checks = {
    "new_segments_vector": sorted(paths.new_segments_vector_dir.glob("??_??_*.gpkg")),
    "bends": sorted(paths.bends_dir.glob("??_??_*.parquet")),
    "single_values": sorted(paths.single_values_dir.glob("??_??_*")),
    "reach_averaged": sorted(paths.reach_averaged_dir.glob("??_??_*")),
    "single_smoothed": sorted(paths.single_smoothed_dir.glob("??_*_smoothed.nc")),
}

foreign = []
for label, files in checks.items():
    for path in files:
        prefix = path.name.split("_")[0]
        if len(prefix) == 2 and prefix != target_continent:
            foreign.append(f"{label}: {path}")

if foreign:
    print(
        f"Refusing to run continent test for '{target_continent}' because the configured results_root "
        "already contains outputs for other continents."
    )
    print("Use a clean test config/results folder or remove the mixed outputs first.")
    for item in foreign:
        print(item)
    raise SystemExit(2)
PY

echo
echo "[1/8] Step 1.5 FABDEM index"
python -m pipeline.build_fabdem_index --config "$CONFIG_PATH"

echo
echo "[2/8] Step 1 segmentation for ${CONTINENT}"
python -m pipeline.segment_reaches "$CONTINENT" --workers 1 --config "$CONFIG_PATH"

echo
echo "[3/8] Step 2 canonical bend-table build for ${CONTINENT}"
python -m pipeline.build_step2_results \
  "$CONTINENT" "$PROCESSORS" \
  --config "$CONFIG_PATH"

echo
echo "[4/8] Step 4 confinement-factor table for ${CONTINENT}"
python -m pipeline.build_step4_confinement_factor \
  "$CONTINENT" \
  --config "$CONFIG_PATH"

echo
echo "[5/8] Step 5 confinement outputs for ${CONTINENT}, height factor ${HEIGHT_FACTOR}"
python -m pipeline.build_step5_confinement_outputs \
  "$CONTINENT" "$PROCESSORS" \
  --height-factor "$HEIGHT_FACTOR" \
  --config "$CONFIG_PATH"

echo
echo "[6/8] Step 6 aggregates for height factor ${HEIGHT_FACTOR}"
python -m pipeline.build_step6_aggregates \
  --height-factor "$HEIGHT_FACTOR" \
  --config "$CONFIG_PATH"

echo
echo "[7/8] Step 7 spatial smoothing for ${CONTINENT}, height factor ${HEIGHT_FACTOR}"
python -m pipeline.build_step7_spatial_smoothing \
  --height-factor "$HEIGHT_FACTOR" \
  --workers 1 \
  --continents "$CONTINENT" \
  --config "$CONFIG_PATH"

echo
echo "[8/8] Step 8 clustering smoke test for height factor ${HEIGHT_FACTOR}"
python -m pipeline.build_step8_confinement_clustering \
  --height-factors "$HEIGHT_FACTOR" \
  --clusters 2 3 \
  --sample-size 10 \
  --workers 1 \
  --config "$CONFIG_PATH"

echo
echo "Continent test finished for ${CONTINENT}."
echo "For a full final-workflow run after building Step 5 and Step 6 for height factors 2 3 4:"
echo "  python -m pipeline.build_step7_8_confinement_workflow --height-factors 2 3 4 --continents ${CONTINENT} --smoothing-workers 1 --clustering-workers 1 --clusters 2 3 --sample-size 10 --config $CONFIG_PATH"
