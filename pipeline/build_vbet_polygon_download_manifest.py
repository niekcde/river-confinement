"""Export selected audit DGO identifiers for targeted Riverscapes project retrieval."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd


AUDIT_GPKG = Path("results/vbet_validation/crosswalk_audit/crosswalk_audit.gpkg")
OUTPUT = Path("results/vbet_validation/crosswalk_audit/selected_dgo_download_manifest.csv")


def main() -> None:
    selected = gpd.read_file(AUDIT_GPKG, layer="selected_dgo_centroids").to_crs(4326)
    columns = [
        "bendID", "reach_id", "audit_stratum", "dgoid", "level_path", "fcode",
        "seg_distance", "match_distance_m", "sword_bend_width_m", "vbet_channel_width",
        "candidate_channel_width_to_sword", "integrated_width",
    ]
    manifest = selected[columns].copy()
    manifest["longitude"] = selected.geometry.x
    manifest["latitude"] = selected.geometry.y
    manifest["project_lookup"] = "Use the DGO ID and coordinates to identify the source VBET/RME project in Riverscapes Data Exchange."
    manifest = manifest.sort_values(["audit_stratum", "bendID", "dgoid"])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(OUTPUT, index=False)
    print(f"Wrote {len(manifest)} selected DGO records for {manifest.bendID.nunique()} audit bends: {OUTPUT}")


if __name__ == "__main__":
    main()
