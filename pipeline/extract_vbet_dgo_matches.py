"""Extract only VBET/RME DGOs near manuscript bend geometries and match them.

Riverscapes' Custom Metrics Dataset may contain tens of millions of DGO
centroid points. This utility queries its GeoPackage spatial index tile by tile,
then performs an exact nearest match in a metre-based CRS. It avoids loading the
complete RME export into memory.
"""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import from_wkb
from shapely.geometry import box


DEFAULT_CRS = "EPSG:5070"


def _gpkg_wkb_to_geometry(value: bytes):
    """Strip a GeoPackage binary header and decode its embedded WKB geometry."""
    if value[:2] != b"GP":
        raise ValueError("Expected a GeoPackage binary geometry.")
    flags = value[3]
    envelope_code = (flags >> 1) & 0b111
    envelope_sizes = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}
    if envelope_code not in envelope_sizes:
        raise ValueError(f"Unsupported GeoPackage envelope code: {envelope_code}")
    return from_wkb(value[8 + envelope_sizes[envelope_code] :])


def _active_tiles(bends: gpd.GeoDataFrame, tile_degrees: float, degree_pad: float):
    west, south, east, north = bends.total_bounds
    x_values = np.arange(np.floor(west / tile_degrees) * tile_degrees, east + tile_degrees, tile_degrees)
    y_values = np.arange(np.floor(south / tile_degrees) * tile_degrees, north + tile_degrees, tile_degrees)
    spatial_index = bends.sindex
    for x_value in x_values:
        for y_value in y_values:
            tile = box(x_value, y_value, x_value + tile_degrees, y_value + tile_degrees)
            expanded = tile.buffer(degree_pad)
            if len(spatial_index.query(expanded, predicate="intersects")):
                yield expanded.bounds


def _fetch_dgos(
    connection: sqlite3.Connection, bounds: tuple[float, float, float, float]
) -> list[tuple]:
    west, south, east, north = bounds
    return connection.execute(
        """
        SELECT d.dgoid, d.geom, g.integrated_width, g.channel_width,
               g.confinement_ratio, g.constriction_ratio, d.level_path,
               d.seg_distance, d.centerline_length, d.segment_area
        FROM rtree_dgos_geom AS r
        JOIN dgos AS d ON d.dgoid = r.id
        JOIN dgo_geomorph AS g ON g.dgoid = d.dgoid
        WHERE r.minx <= ? AND r.maxx >= ? AND r.miny <= ? AND r.maxy >= ?
        """,
        (east, west, north, south),
    ).fetchall()


def match_dgos_to_bends(
    dgo_gpkg: str | Path,
    bends_path: str | Path,
    output_path: str | Path,
    *,
    bend_id: str,
    search_distance_m: float = 250.0,
    tile_degrees: float = 0.25,
    projected_crs: str = DEFAULT_CRS,
    bends_layer: str | None = None,
) -> Path:
    if search_distance_m <= 0 or tile_degrees <= 0:
        raise ValueError("Search distance and tile size must both be positive.")
    connection = sqlite3.connect(dgo_gpkg)
    dgo_bounds = connection.execute(
        "SELECT min_x, min_y, max_x, max_y FROM gpkg_contents WHERE table_name = 'dgos'"
    ).fetchone()
    if dgo_bounds is None:
        raise ValueError("The supplied GeoPackage has no registered 'dgos' layer.")
    if any(value is None for value in dgo_bounds):
        dgo_bounds = connection.execute(
            "SELECT min(minx), min(miny), max(maxx), max(maxy) FROM rtree_dgos_geom"
        ).fetchone()
    # Read only the manuscript bends within the RME-export envelope.  Reading a
    # global bend GeoPackage in full is unnecessary and can dominate run time.
    bends = gpd.read_file(bends_path, layer=bends_layer, bbox=dgo_bounds)
    if bends.crs is None:
        raise ValueError("Manuscript bend layer has no CRS.")
    if bend_id not in bends:
        raise KeyError(f"Manuscript bend layer has no {bend_id!r} column.")
    bends = bends.loc[bends.geometry.notna() & ~bends.geometry.is_empty].copy()
    bends_wgs84 = bends.to_crs("EPSG:4326")

    # At the northern edge of CONUS, one degree of longitude is about 74 km.
    # This conservative pad ensures a 250 m projected-distance match is never
    # missed while still making the GeoPackage R-tree do the heavy filtering.
    degree_pad = search_distance_m / 70_000.0
    candidates: dict[int, tuple] = {}
    for tile_bounds in _active_tiles(bends_wgs84, tile_degrees, degree_pad):
        for row in _fetch_dgos(connection, tile_bounds):
            candidates[row[0]] = row
    connection.close()
    if not candidates:
        raise RuntimeError("No RME DGOs were found near the supplied bend geometries.")

    columns = [
        "dgoid", "geom", "integrated_width", "vbet_channel_width",
        "vbet_confinement_ratio", "vbet_constriction_ratio", "level_path",
        "seg_distance", "centerline_length", "segment_area",
    ]
    values = list(candidates.values())
    dgos = gpd.GeoDataFrame(
        pd.DataFrame(values, columns=columns).drop(columns="geom"),
        geometry=[_gpkg_wkb_to_geometry(row[1]) for row in values],
        crs="EPSG:4326",
    ).to_crs(projected_crs)
    bends_projected = bends.to_crs(projected_crs)
    matches = gpd.sjoin_nearest(
        dgos,
        bends_projected,
        how="left",
        max_distance=search_distance_m,
        distance_col="match_distance_m",
        lsuffix="vbet",
        rsuffix="bend",
    )
    tie_counts = matches.groupby("dgoid").size().rename("nearest_tie_count")
    matches = matches.join(tie_counts, on="dgoid")
    matches["match_status"] = "matched"
    matches.loc[matches[bend_id].isna(), "match_status"] = "no_bend_within_threshold"
    matches.loc[
        matches["nearest_tie_count"].gt(1) & matches[bend_id].notna(), "match_status"
    ] = "ambiguous_nearest_tie"
    matches = matches.loc[matches["match_status"].eq("matched")].copy()
    matches["vbet_width_to_vbet_channel_ratio"] = (
        matches["integrated_width"] / matches["vbet_channel_width"]
    )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    matches.to_file(output, driver="GPKG")
    matches.drop(columns="geometry").to_parquet(output.with_suffix(".parquet"), index=False)
    return output


def main_cli(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Match VBET/RME DGO centroids to manuscript bends.")
    parser.add_argument("--dgos", required=True, help="Riverscapes riverscape_metrics.gpkg path")
    parser.add_argument("--bends", required=True, help="Manuscript bend GeoPackage path")
    parser.add_argument("--output", required=True)
    parser.add_argument("--bend-id", default="bendID")
    parser.add_argument("--search-distance-m", type=float, default=250.0)
    parser.add_argument("--tile-degrees", type=float, default=0.25)
    parser.add_argument("--projected-crs", default=DEFAULT_CRS)
    parser.add_argument("--bends-layer")
    args = parser.parse_args(argv)
    output = match_dgos_to_bends(
        args.dgos, args.bends, args.output, bend_id=args.bend_id,
        search_distance_m=args.search_distance_m, tile_degrees=args.tile_degrees,
        projected_crs=args.projected_crs, bends_layer=args.bends_layer,
    )
    print(output)


if __name__ == "__main__":
    main_cli()
