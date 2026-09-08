"""Resolve coarse VBET-to-SWORD candidate matches at the VBET level-path scale.

The 150 m centroid match is a candidate generator only. This script groups
candidate DGOs by VBET level_path and selects the path whose VBET channel width
is most compatible with the manuscript SWORD bend width. NHD FCode is retained
as a tie-breaker and quality flag, not a universal exclusion criterion.
"""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd


def args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matches", type=Path, default=Path("results/vbet_validation/pnw_dgo_bend_matches.parquet"))
    parser.add_argument("--bend-values", type=Path, default=Path("results/vbet_validation/pnw_bend_height_factor_sensitivity.parquet"))
    parser.add_argument(
        "--vbet-gpkg",
        type=Path,
        default=Path("data/vbet/report_Custom_Riverscapes_Metrics_Dataset_Custom_Riverscapes_Metrics_Dataset_Report_5d9b73d4-43ee-4ef6-91b8-1731ada38c6e/outputs/riverscape_metrics.gpkg"),
    )
    parser.add_argument("--output", type=Path, default=Path("results/vbet_validation/pnw_refined_vbet_level_path_matches.parquet"))
    parser.add_argument("--match-distance-m", type=float, default=150)
    return parser.parse_args()


def fetch_fcodes(gpkg: Path, dgoids: pd.Series) -> pd.DataFrame:
    """Fetch FCode only for supplied DGOs, without scanning the 8.9 GB layer."""
    ids = pd.DataFrame({"dgoid": pd.unique(dgoids).astype("int64")})
    with sqlite3.connect(gpkg) as con:
        con.execute("CREATE TEMP TABLE selected_dgos (dgoid INTEGER PRIMARY KEY)")
        con.executemany("INSERT INTO selected_dgos VALUES (?)", ((int(value),) for value in ids.dgoid))
        return pd.read_sql_query(
            "SELECT d.dgoid, d.fcode FROM dgos AS d INNER JOIN selected_dgos AS s ON d.dgoid = s.dgoid",
            con,
        )


def fcode_priority(fcode: object) -> int:
    """Use intermittent status only to break otherwise equal width matches."""
    return 1 if fcode == 46003 else 0


def main() -> None:
    options = args()
    for path in (options.matches, options.bend_values, options.vbet_gpkg):
        if not path.exists():
            raise FileNotFoundError(path)

    matches = pd.read_parquet(options.matches)
    matches = matches.loc[matches.match_distance_m <= options.match_distance_m].copy()
    bends = pd.read_parquet(options.bend_values)[["bendID", "bendWidths"]].drop_duplicates("bendID")
    matches = matches.merge(bends, on="bendID", how="inner", validate="many_to_one")
    matches = matches.merge(fetch_fcodes(options.vbet_gpkg, matches.dgoid), on="dgoid", how="left", validate="many_to_one")

    grouped = (
        matches.groupby(["bendID", "level_path"], as_index=False)
        .agg(
            n_dgos_on_level_path=("dgoid", "size"),
            dgoids=("dgoid", lambda values: ",".join(map(str, sorted(values.unique())))),
            fcodes=("fcode", lambda values: ",".join(str(int(value)) for value in sorted(values.dropna().unique()))),
            vbet_integrated_width_m=("integrated_width", "median"),
            vbet_channel_width_m=("vbet_channel_width", "median"),
            median_match_distance_m=("match_distance_m", "median"),
            min_match_distance_m=("match_distance_m", "min"),
            bendWidths=("bendWidths", "first"),
        )
    )
    grouped["channel_width_ratio_to_sword"] = grouped.vbet_channel_width_m / grouped.bendWidths
    positive_width_ratio = grouped.channel_width_ratio_to_sword.where(grouped.channel_width_ratio_to_sword > 0)
    grouped["channel_width_log_difference"] = np.abs(np.log10(positive_width_ratio)).fillna(np.inf)
    grouped["fcode_priority"] = grouped.fcodes.map(lambda value: fcode_priority(int(value)) if value.isdigit() else 0)
    grouped["n_candidate_level_paths"] = grouped.groupby("bendID").level_path.transform("size")

    # FCode is intentionally after width compatibility in the sorting key.
    grouped = grouped.sort_values(
        ["bendID", "channel_width_log_difference", "fcode_priority", "median_match_distance_m"],
        kind="stable",
    )
    grouped["candidate_rank"] = grouped.groupby("bendID").cumcount() + 1
    grouped["next_best_width_score"] = grouped.groupby("bendID").channel_width_log_difference.shift(-1)
    grouped["width_score_separation"] = grouped.next_best_width_score - grouped.channel_width_log_difference
    selected = grouped.loc[grouped.candidate_rank == 1].copy()
    selected["channel_width_compatibility"] = pd.cut(
        selected.channel_width_log_difference,
        bins=[-np.inf, np.log10(3), 1, np.inf],
        labels=["within factor 3", "within factor 10", "outside factor 10"],
    )
    selected["fcode_interpretation"] = np.where(
        selected.fcodes.eq("46003"), "intermittent-stream candidate", "not solely intermittent-stream candidate"
    )
    selected["crosswalk_tier"] = np.select(
        [
            selected.channel_width_ratio_to_sword.between(0.5, 2.0),
            selected.channel_width_ratio_to_sword.between(0.2, 5.0),
        ],
        ["high confidence: within factor 2", "sensitivity only: factor 2 to 5"],
        default="reject: outside factor 5",
    )
    selected["accepted_high_confidence"] = selected.channel_width_ratio_to_sword.between(0.5, 2.0)
    selected["accepted_relaxed_sensitivity"] = selected.channel_width_ratio_to_sword.between(0.2, 5.0)

    options.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_parquet(options.output, index=False)
    print(f"Wrote {len(selected):,} selected level-path candidates to {options.output}")
    print(selected.channel_width_compatibility.value_counts(dropna=False).to_string())
    print(selected.crosswalk_tier.value_counts(dropna=False).to_string())
    print("Bends with >1 candidate level path:", int((selected.n_candidate_level_paths > 1).sum()))


if __name__ == "__main__":
    main()
