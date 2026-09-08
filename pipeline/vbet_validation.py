"""Spatial crosswalk and width comparison for SWORD and VBET/RME IGO exports.

This module intentionally does not download Riverscapes projects.  Downloading is
best handled by the Riverscapes Reports interface (for an AOI) or its API; this
tool starts with the resulting IGO point export.  That keeps the validation
workflow independent of multi-terabyte VBET polygon archives.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd


DEFAULT_METRE_CRS = "EPSG:5070"  # CONUS Albers Equal Area


def _read_vector(path: str | Path, layer: str | None = None) -> gpd.GeoDataFrame:
    frame = gpd.read_file(path, layer=layer)
    if frame.crs is None:
        raise ValueError(f"{path} has no CRS; assign its CRS before spatial matching.")
    frame = frame.loc[frame.geometry.notna() & ~frame.geometry.is_empty].copy()
    return frame


def _project(frame: gpd.GeoDataFrame, crs: str) -> gpd.GeoDataFrame:
    return frame.to_crs(crs) if str(frame.crs) != crs else frame.copy()


def build_candidate_universe(
    sword_path: str | Path,
    output_path: str | Path,
    *,
    reach_id: str,
    metric: str,
    min_metric: float | None = None,
    max_metric: float | None = None,
    layer: str | None = None,
) -> Path:
    """Write the valid SWORD records that are eligible for a VBET match.

    A VBET export, rather than an assumed national footprint, determines actual
    coverage.  This function therefore applies only validity filters that can
    be evaluated from SWORD/manuscript data.
    """
    sword = _read_vector(sword_path, layer)
    for column in (reach_id, metric):
        if column not in sword:
            raise KeyError(f"SWORD input is missing required column {column!r}.")
    metric_values = pd.to_numeric(sword[metric], errors="coerce")
    valid = sword[reach_id].notna() & np.isfinite(metric_values)
    if min_metric is not None:
        valid &= metric_values >= min_metric
    if max_metric is not None:
        valid &= metric_values <= max_metric
    candidate = sword.loc[valid].copy()
    candidate[metric] = metric_values.loc[valid]
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    candidate.to_file(output, driver="GPKG")
    return output


def match_igos_to_sword(
    sword_path: str | Path,
    igo_path: str | Path,
    output_path: str | Path,
    *,
    reach_id: str,
    igo_id: str,
    max_distance_m: float = 250.0,
    projected_crs: str = DEFAULT_METRE_CRS,
    sword_layer: str | None = None,
    igo_layer: str | None = None,
) -> Path:
    """Match each VBET IGO point to its nearest eligible SWORD reach geometry.

    The output retains the matching distance and a tie flag.  Records at or
    beyond ``max_distance_m`` and exact nearest-distance ties are excluded from
    the primary match set; retain the input/export for visual QA of exclusions.
    """
    if max_distance_m <= 0:
        raise ValueError("max_distance_m must be positive.")
    sword = _project(_read_vector(sword_path, sword_layer), projected_crs)
    igos = _project(_read_vector(igo_path, igo_layer), projected_crs)
    for frame, column, label in ((sword, reach_id, "SWORD"), (igos, igo_id, "IGO")):
        if column not in frame:
            raise KeyError(f"{label} input is missing ID column {column!r}.")
        if frame[column].isna().any():
            raise ValueError(f"{label} ID column {column!r} contains null values.")
        if not frame[column].is_unique:
            raise ValueError(f"{label} ID column {column!r} must be unique for matching.")

    # Keep SWORD attributes (especially the channel/bend-width denominator) in
    # the crosswalk so aggregation does not need a second, potentially lossy
    # spatial join.  If an IGO export happens to use a duplicate column name,
    # GeoPandas applies its documented suffixes; users should then pass that
    # resulting field name to ``aggregate --denominator``.
    igos = igos.copy()
    joined = gpd.sjoin_nearest(
        igos,
        sword,
        how="left",
        max_distance=max_distance_m,
        distance_col="match_distance_m",
    )
    counts = joined.groupby(igo_id, dropna=False).size().rename("nearest_tie_count")
    joined = joined.join(counts, on=igo_id)
    joined["match_status"] = "matched"
    joined.loc[joined[reach_id].isna(), "match_status"] = "no_reach_within_threshold"
    joined.loc[
        joined["nearest_tie_count"].gt(1) & joined[reach_id].notna(), "match_status"
    ] = "ambiguous_nearest_tie"
    matched = joined.loc[joined["match_status"].eq("matched")].copy()
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    matched.to_file(output, driver="GPKG")
    summary = (
        joined.groupby("match_status", dropna=False)
        .size()
        .rename("igo_count")
        .reset_index()
    )
    summary.to_csv(output.with_suffix(".match_summary.csv"), index=False)
    return output


def aggregate_vbet_widths(
    matches_path: str | Path,
    output_path: str | Path,
    *,
    reach_id: str,
    width: str,
    denominator: str,
    min_igos: int = 3,
    layer: str | None = None,
) -> Path:
    """Aggregate matched IGO widths and construct the VBET/SWORD width ratio.

    The primary aggregation is the median IGO integrated width per SWORD reach,
    which resists a single anomalously broad valley-bottom segment.  Length-
    weighted means can be added later if the report export contains a clearly
    documented DGO/IGO centreline-length field.
    """
    matches = _read_vector(matches_path, layer)
    for column in (reach_id, width, denominator):
        if column not in matches:
            raise KeyError(f"Matched input is missing required column {column!r}.")
    matches[width] = pd.to_numeric(matches[width], errors="coerce")
    matches[denominator] = pd.to_numeric(matches[denominator], errors="coerce")
    valid = np.isfinite(matches[width]) & np.isfinite(matches[denominator]) & (matches[denominator] > 0)
    grouped = matches.loc[valid].groupby(reach_id, as_index=False).agg(
        vbet_igo_count=(width, "count"),
        vbet_width_m=(width, "median"),
        vbet_width_iqr_m=(width, lambda values: values.quantile(0.75) - values.quantile(0.25)),
        sword_denominator_m=(denominator, "median"),
        median_match_distance_m=("match_distance_m", "median"),
    )
    grouped = grouped.loc[grouped["vbet_igo_count"] >= min_igos].copy()
    grouped["vbet_to_sword_width_ratio"] = grouped["vbet_width_m"] / grouped["sword_denominator_m"]
    grouped["log10_vbet_to_sword_width_ratio"] = np.log10(grouped["vbet_to_sword_width_ratio"])
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    grouped.to_parquet(output, index=False)
    return output


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build a spatial SWORD–VBET width comparison.")
    commands = parser.add_subparsers(dest="command", required=True)
    candidate = commands.add_parser("candidates", help="Filter valid SWORD records for matching.")
    candidate.add_argument("--sword", required=True)
    candidate.add_argument("--output", required=True)
    candidate.add_argument("--reach-id", required=True)
    candidate.add_argument("--metric", required=True, help="Existing manuscript confinement metric column.")
    candidate.add_argument("--min-metric", type=float)
    candidate.add_argument("--max-metric", type=float)
    candidate.add_argument("--layer")
    match = commands.add_parser("match", help="Spatially match VBET IGO points to SWORD reaches.")
    match.add_argument("--sword", required=True)
    match.add_argument("--igos", required=True)
    match.add_argument("--output", required=True)
    match.add_argument("--reach-id", required=True)
    match.add_argument("--igo-id", required=True)
    match.add_argument("--max-distance-m", type=float, default=250)
    match.add_argument("--projected-crs", default=DEFAULT_METRE_CRS)
    match.add_argument("--sword-layer")
    match.add_argument("--igo-layer")
    aggregate = commands.add_parser("aggregate", help="Aggregate IGO widths by matched SWORD reach.")
    aggregate.add_argument("--matches", required=True)
    aggregate.add_argument("--output", required=True)
    aggregate.add_argument("--reach-id", required=True)
    aggregate.add_argument("--width", required=True, help="VBET integrated-width column from the export.")
    aggregate.add_argument("--denominator", required=True, help="SWORD channel/bend-width column.")
    aggregate.add_argument("--min-igos", type=int, default=3)
    aggregate.add_argument("--layer")
    return parser


def main_cli(argv: list[str] | None = None) -> None:
    args = _parser().parse_args(argv)
    if args.command == "candidates":
        output = build_candidate_universe(args.sword, args.output, reach_id=args.reach_id, metric=args.metric, min_metric=args.min_metric, max_metric=args.max_metric, layer=args.layer)
    elif args.command == "match":
        output = match_igos_to_sword(args.sword, args.igos, args.output, reach_id=args.reach_id, igo_id=args.igo_id, max_distance_m=args.max_distance_m, projected_crs=args.projected_crs, sword_layer=args.sword_layer, igo_layer=args.igo_layer)
    else:
        output = aggregate_vbet_widths(args.matches, args.output, reach_id=args.reach_id, width=args.width, denominator=args.denominator, min_igos=args.min_igos, layer=args.layer)
    print(output)


if __name__ == "__main__":
    main_cli()
