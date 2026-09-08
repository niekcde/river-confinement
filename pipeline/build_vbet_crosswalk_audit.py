"""Create a stratified, QGIS-ready audit sample for the VBET-SWORD crosswalk.

The Riverscapes report export contains DGO centroids rather than DGO polygons.
This package therefore supports river-identity review with candidate centroids,
SWORD bend locations, and full candidate attributes; it does not claim polygon
overlap verification.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refined", type=Path, default=Path("results/vbet_validation/pnw_refined_vbet_level_path_matches.parquet"))
    parser.add_argument("--factor-results", type=Path, default=Path("results/vbet_validation/pnw_bend_height_factor_sensitivity.parquet"))
    parser.add_argument("--coarse-matches", type=Path, default=Path("results/vbet_validation/pnw_dgo_bend_matches.parquet"))
    parser.add_argument(
        "--vbet-gpkg",
        type=Path,
        default=Path("data/vbet/report_Custom_Riverscapes_Metrics_Dataset_Custom_Riverscapes_Metrics_Dataset_Report_5d9b73d4-43ee-4ef6-91b8-1731ada38c6e/outputs/riverscape_metrics.gpkg"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/vbet_validation/crosswalk_audit"))
    parser.add_argument("--per-stratum", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--anchor-bend-ids", nargs="*", default=["12005_3", "12132_9"])
    return parser.parse_args()


def take(frame: pd.DataFrame, n: int, seed: int, *, random: bool = True) -> pd.DataFrame:
    if random:
        return frame.sample(min(n, len(frame)), random_state=seed)
    return frame.head(n)


def read_dgo_centroids(gpkg: Path, candidates: pd.DataFrame) -> gpd.GeoDataFrame:
    """Read only candidate centroid features, then restore DGO IDs via stable keys."""
    ids = candidates.dgoid.drop_duplicates().astype(int).tolist()
    batches = []
    for start in range(0, len(ids), 800):
        batch = ids[start : start + 800]
        where = "dgoid IN ({})".format(",".join(map(str, batch)))
        batches.append(gpd.read_file(gpkg, layer="dgos", where=where))
    dgos = pd.concat(batches, ignore_index=True)
    dgos = gpd.GeoDataFrame(dgos, geometry="geometry", crs=batches[0].crs)
    # dgoid is the GPKG feature id and not exposed by GeoPandas. level_path +
    # seg_distance is sufficient for this audit export and joins back to IDs.
    lookup = candidates[["dgoid", "level_path", "seg_distance"]].drop_duplicates()
    return dgos.merge(lookup, on=["level_path", "seg_distance"], how="inner", validate="one_to_one")


def main() -> None:
    args = parse_args()
    for path in (args.refined, args.factor_results, args.coarse_matches, args.vbet_gpkg):
        if not path.exists():
            raise FileNotFoundError(path)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    refined = pd.read_parquet(args.refined)
    factors = pd.read_parquet(args.factor_results)
    audit = refined.merge(
        factors[["bendID", "n_vbet_dgos", "paper_ratio_h2", "vbet_ratio"]],
        on="bendID", how="left", validate="one_to_one",
    )
    audit["refined_vbet_ratio"] = audit.vbet_integrated_width_m / audit.bendWidths
    audit["log10_refined_vbet_over_paper"] = np.log10(audit.refined_vbet_ratio / audit.paper_ratio_h2)
    audit = audit.loc[(audit.n_vbet_dgos >= 3) & np.isfinite(audit.log10_refined_vbet_over_paper)].copy()

    # Keep categories disjoint, retaining two known illustrative cases.
    pieces = []
    anchors = audit.loc[audit.bendID.isin(args.anchor_bend_ids)].copy()
    anchors["audit_stratum"] = "named inspection case"
    pieces.append(anchors)
    used = set(anchors.bendID)
    remaining = audit.loc[~audit.bendID.isin(used)].sort_values("log10_refined_vbet_over_paper")
    extreme = take(remaining, args.per_stratum, args.seed, random=False).assign(audit_stratum="largest remaining paper-over-VBET discrepancy")
    pieces.append(extreme)
    used.update(extreme.bendID)
    strata = [
        ("outside factor-10 channel-width compatibility", audit.channel_width_compatibility.eq("outside factor 10")),
        ("multiple VBET level-path candidates", audit.n_candidate_level_paths.gt(1)),
        ("single, width-compatible candidate", audit.n_candidate_level_paths.eq(1) & audit.channel_width_compatibility.eq("within factor 3")),
    ]
    for i, (name, condition) in enumerate(strata, start=1):
        pool = audit.loc[condition & ~audit.bendID.isin(used)]
        chosen = take(pool, args.per_stratum, args.seed + i).assign(audit_stratum=name)
        pieces.append(chosen)
        used.update(chosen.bendID)
    sample = pd.concat(pieces, ignore_index=True).sort_values(["audit_stratum", "bendID"])
    sample.to_csv(args.output_dir / "crosswalk_audit_sample.csv", index=False)

    coarse = pd.read_parquet(args.coarse_matches)
    candidates = coarse.loc[(coarse.match_distance_m <= 150) & coarse.bendID.isin(sample.bendID)].copy()
    selected_ids = sample[["bendID", "level_path"]].rename(columns={"level_path": "selected_level_path"})
    candidates = candidates.merge(selected_ids, on="bendID", how="left", validate="many_to_one")
    candidates["is_selected_level_path"] = candidates.level_path.eq(candidates.selected_level_path)
    dgo_centroids = read_dgo_centroids(args.vbet_gpkg, candidates)
    candidate_points = candidates.drop(columns="geometry", errors="ignore").merge(
        dgo_centroids[["dgoid", "fcode", "geometry"]], on="dgoid", how="left", validate="many_to_one"
    )
    candidate_points = gpd.GeoDataFrame(candidate_points, geometry="geometry", crs=dgo_centroids.crs)
    audit_attributes = sample[[
        "bendID", "audit_stratum", "bendWidths", "vbet_channel_width_m",
        "channel_width_ratio_to_sword", "channel_width_compatibility",
        "n_candidate_level_paths", "width_score_separation",
    ]].rename(columns={
        "bendWidths": "sword_bend_width_m",
        "vbet_channel_width_m": "selected_vbet_channel_width_m",
    })
    candidate_points = candidate_points.merge(audit_attributes, on="bendID", how="left", validate="many_to_one")
    candidate_points["candidate_channel_width_to_sword"] = (
        candidate_points["vbet_channel_width"] / candidate_points["sword_bend_width_m"]
    )

    bend_points = gpd.read_file(args.coarse_matches.with_suffix(".gpkg"))
    bend_points = bend_points.loc[bend_points.bendID.isin(sample.bendID)].drop_duplicates("bendID")
    bend_points = bend_points.merge(sample, on="bendID", how="left", validate="one_to_one")
    output_gpkg = args.output_dir / "crosswalk_audit.gpkg"
    bend_points.to_file(output_gpkg, layer="audit_bends", driver="GPKG")
    candidate_points.to_file(output_gpkg, layer="candidate_dgo_centroids", driver="GPKG")
    candidate_points.loc[candidate_points.is_selected_level_path].to_file(output_gpkg, layer="selected_dgo_centroids", driver="GPKG")
    print(f"Wrote {len(sample)} audit bends and {len(candidate_points)} candidate DGO centroids to {output_gpkg}")


if __name__ == "__main__":
    main()
