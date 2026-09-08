#!/usr/bin/env python3
"""Geometric crosswalk audit between one downloaded VBET project and SWORD bends.

This is deliberately an *identity/plausibility* diagnostic, not another valley
width estimate.  It asks whether several adjacent SWORD bend centre lines follow
the VBET ``vbet_full`` valley-bottom corridor selected by the DGO/level-path
crosswalk.

Example
-------
python pipeline/check_vbet_sword_polygon_overlap.py \
  --dgoid 9753974 \
  --vbet-gpkg results/vbet_validation/crosswalk_audit/audits/vbet_9753974.gpkg

The default source SWORD bends are the submitted manuscript product on the PhD
volume.  Override ``--sword-bends`` if it is mounted elsewhere.
"""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pyogrio
from shapely.ops import unary_union


DEFAULT_AUDIT = Path("results/vbet_validation/crosswalk_audit/crosswalk_audit.gpkg")
DEFAULT_SWORD = Path("/Volumes/PhD/confinement/results/submission/v1/global_clusters_bend.gpkg")
DEFAULT_OUTPUT = Path("results/vbet_validation/crosswalk_audit/polygon_overlap_checks.csv")
METRIC_CRS = "EPSG:5070"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dgoid", required=True, type=int, help="VBET DGO ID in candidate_dgo_centroids")
    parser.add_argument("--vbet-gpkg", required=True, type=Path, help="Downloaded VBET output GeoPackage")
    parser.add_argument("--audit-gpkg", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--audit-path", choices=("coarse", "refined"), default="coarse",
                        help=("For a DGO read from audit_bends, test its original coarse path "
                              "(default) or the width-refined selected path"))
    parser.add_argument("--sword-bends", type=Path, default=DEFAULT_SWORD)
    parser.add_argument("--neighbour-bends", type=int, default=2,
                        help="Number of SWORD bends on either side of target bend (default: 2)")
    parser.add_argument("--corridor-buffer-m", type=float, default=10.0,
                        help="Buffer used only for alignment sensitivity metric (default: 10 m)")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    candidates = gpd.read_file(args.audit_gpkg, layer="candidate_dgo_centroids")
    hit = candidates.loc[candidates["dgoid"].eq(args.dgoid)]
    if len(hit) == 1:
        hit = hit.iloc[0]
        level_path = str(hit["level_path"])
        target_rank = int(hit["bendRank"])
        reach_id = round(float(hit["reach_id"]))
        target_bend_id = hit["bendID"]
        path_source = "candidate_dgo_centroids"
    else:
        # audit_bends preserves the original coarse nearest-DGO record in its
        # dgoid/level_path_x columns, even where that point is not in the later
        # candidate-centroid layer.  This is the ID users see when inspecting the
        # audit-bends layer in QGIS.
        audit_bends = gpd.read_file(args.audit_gpkg, layer="audit_bends")
        audit_hit = audit_bends.loc[audit_bends["dgoid"].eq(args.dgoid)]
        if len(audit_hit) == 1:
            audit_hit = audit_hit.iloc[0]
            level_column = "level_path_x" if args.audit_path == "coarse" else "level_path_y"
            level_path = str(audit_hit[level_column])
            target_rank = int(audit_hit["bendRank"])
            reach_id = round(float(audit_hit["reach_id"]))
            target_bend_id = audit_hit["bendID"]
            path_source = f"audit_bends ({args.audit_path})"
        else:
            # A downloaded VBET GeoPackage does not retain the Report DGO ID.
            # List the Report candidates spatially covered by this download.
            with sqlite3.connect(args.vbet_gpkg) as connection:
                bounds = connection.execute(
                    "SELECT min_x, min_y, max_x, max_y FROM gpkg_contents "
                    "WHERE table_name = 'vbet_igos'"
                ).fetchone()
            if bounds is None:
                raise ValueError(f"Could not find the vbet_igos extent in {args.vbet_gpkg}")
            minx, miny, maxx, maxy = bounds
            in_download = candidates.to_crs(4326).cx[minx:maxx, miny:maxy]
            columns = ["dgoid", "bendID", "level_path", "seg_distance", "is_selected_level_path"]
            available = in_download[columns].sort_values(["bendID", "dgoid"]).to_string(index=False)
            raise ValueError(
                f"DGO {args.dgoid} is not present in either audit layer.\n\n"
                f"Candidate DGO records spatially covered by {args.vbet_gpkg.name}:\n{available}"
            )

    # Only the selected VBET network path is relevant; nearby tributary polygons
    # must not inflate apparent overlap.
    vbet = pyogrio.read_dataframe(
        args.vbet_gpkg,
        layer="vbet_full",
        where=f"level_path = '{level_path}'",
    )
    if vbet.empty:
        raise ValueError(f"No vbet_full polygons for level_path {level_path}")

    # Read only the small geographic neighbourhood covered by this VBET project,
    # then retain the matched SWORD reach and a local multi-bend window.
    vbet_ll = vbet.to_crs(4326)
    minx, miny, maxx, maxy = vbet_ll.total_bounds
    pad = 0.002  # protects against tiny cross-dataset positional offsets
    sword = pyogrio.read_dataframe(
        args.sword_bends,
        layer="global_clusters_bend",
        bbox=(minx - pad, miny - pad, maxx + pad, maxy + pad),
        columns=["bendID", "bendRank", "reach_id"],
    )
    sword = sword.loc[sword["reach_id"].round().eq(reach_id)].copy()
    if sword.empty:
        raise ValueError(f"No SWORD bends found for reach_id {reach_id} in VBET project extent")
    sword["bendRank"] = sword["bendRank"].astype(int)
    sword = sword.loc[sword["bendRank"].between(target_rank - args.neighbour_bends,
                                                   target_rank + args.neighbour_bends)].sort_values("bendRank")
    if sword.empty:
        raise ValueError("No target-neighbourhood bends remain after rank filter")

    vbet_m = vbet.to_crs(METRIC_CRS)
    sword_m = sword.to_crs(METRIC_CRS)
    valley_bottom = unary_union(vbet_m.geometry.values)

    records = []
    for _, bend in sword_m.iterrows():
        line = bend.geometry
        length_m = line.length
        inside_m = line.intersection(valley_bottom).length
        corridor = line.buffer(args.corridor_buffer_m)
        records.append({
            "dgoid": args.dgoid,
            "level_path": level_path,
            "reach_id": reach_id,
            "target_bendID": target_bend_id,
            "path_source": path_source,
            "bendID": bend["bendID"],
            "bendRank": bend["bendRank"],
            "line_length_m": length_m,
            "line_inside_vbet_m": inside_m,
            "centreline_inside_fraction": inside_m / length_m if length_m else float("nan"),
            "corridor_buffer_m": args.corridor_buffer_m,
            "corridor_inside_fraction": corridor.intersection(valley_bottom).area / corridor.area,
        })

    result = pd.DataFrame(records)
    all_line = unary_union(sword_m.geometry.values)
    total_length = all_line.length
    total_inside = all_line.intersection(valley_bottom).length
    summary = {
        "dgoid": args.dgoid,
        "level_path": level_path,
        "reach_id": reach_id,
        "target_bendID": target_bend_id,
        "path_source": path_source,
        "bendID": "ALL_SELECTED_BENDS",
        "bendRank": pd.NA,
        "line_length_m": total_length,
        "line_inside_vbet_m": total_inside,
        "centreline_inside_fraction": total_inside / total_length,
        "corridor_buffer_m": args.corridor_buffer_m,
        "corridor_inside_fraction": all_line.buffer(args.corridor_buffer_m).intersection(valley_bottom).area /
                                    all_line.buffer(args.corridor_buffer_m).area,
    }
    result = pd.concat([result, pd.DataFrame([summary])], ignore_index=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Preserve results from other downloaded projects while replacing this one.
    if args.output.exists():
        existing = pd.read_csv(args.output)
        existing = existing.loc[existing["dgoid"].ne(args.dgoid)]
        result = pd.concat([existing, result], ignore_index=True)
    result.to_csv(args.output, index=False)
    print(result.to_string(index=False))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
