#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse
import gc
import os
import re
import sqlite3
from itertools import chain
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from osgeo import gdal
from shapely.geometry import LineString

from .get_orthogonals import OrthogonalGeometry, sample_orthogonal_profiles
from .paths import load_project_paths


PROFILE_COLUMNS = [
    "combined_reach_id",
    "centerlineWkt",
    "elevOut",
    "elevInn",
    "distOut",
    "distInn",
    "lineOut",
    "lineInn",
    "LROrthog",
    "lineSlope",
    "bendHeight",
    "sampling_code",
]
ARRAY_STRING_RE = re.compile(r"[()]|list")
ORTHOGONAL_READ_BATCH_SIZE = 20_000


def _current_rss_mb() -> float | None:
    if os.environ.get("STEP2_REPORT_MEMORY") != "1":
        return None
    try:
        import psutil

        return psutil.Process().memory_info().rss / (1024 * 1024)
    except (ImportError, OSError):
        return None


def empty_profile_frame() -> pd.DataFrame:
    return pd.DataFrame(columns=PROFILE_COLUMNS)


def _resolve_conf_factor(orthogonal_path: Path, conf_factor: int | None) -> int:
    if conf_factor is not None:
        return conf_factor

    try:
        return int(orthogonal_path.stem.split("_")[-1])
    except (IndexError, ValueError) as exc:
        raise ValueError(
            f"Unable to infer conf_factor from orthogonal file name: {orthogonal_path.name}"
        ) from exc


def _serialize_profile_array(values):
    if isinstance(values, float) and np.isnan(values):
        return np.nan

    text = np.array2string(values, separator=",")
    return ARRAY_STRING_RE.sub("", text)


def _serialize_line_array(lines):
    if isinstance(lines, float) and np.isnan(lines):
        return np.nan

    return str([geometry.wkt if isinstance(geometry, LineString) else 9999 for geometry in lines])


def _serialize_scalar_array(values):
    if isinstance(values, float) and np.isnan(values):
        return np.nan

    return np.array2string(values, separator=", ")


def _group_value(group, line_type, column):
    values = group.loc[group["line_type"] == line_type, column]
    if values.empty:
        raise ValueError(f"Missing '{line_type}' row for column '{column}'")
    return values.iloc[0]


def _group_geometry(group, line_type):
    geometries = group.loc[group["line_type"] == line_type, "geometry"]
    if geometries.empty:
        raise ValueError(f"Missing '{line_type}' geometry")
    return geometries.iloc[0]


def _reconstruct_combined_reach(group_local: gpd.GeoDataFrame):
    reach_crs = group_local["reach_crs"].dropna().iloc[0]
    centerline = _group_geometry(group_local, "centerline")

    bend_rows = group_local[group_local["bend_index"] >= 0].copy()
    bend_indices = sorted(int(index) for index in bend_rows["bend_index"].dropna().unique())

    bend_lines = np.empty(len(bend_indices), dtype=object)
    line_out = np.empty(len(bend_indices), dtype=object)
    line_inn = np.empty(len(bend_indices), dtype=object)
    bend_widths = np.empty(len(bend_indices))
    left_right = np.empty(len(bend_indices))

    for output_index, bend_index in enumerate(bend_indices):
        bend_group = bend_rows[bend_rows["bend_index"] == bend_index]
        bend_lines[output_index] = _group_geometry(bend_group, "bendline")
        line_out[output_index] = _group_geometry(bend_group, "outer")
        line_inn[output_index] = _group_geometry(bend_group, "inner")
        bend_widths[output_index] = float(_group_value(bend_group, "outer", "bend_width"))
        left_right[output_index] = float(_group_value(bend_group, "outer", "left_right"))

    return reach_crs, centerline, bend_lines, bend_widths, OrthogonalGeometry(line_out, line_inn, left_right)


def _profile_output_path(paths, orthogonal_path: Path, output_path: str | None = None) -> Path:
    if output_path is not None:
        output = Path(output_path).expanduser()
        if not output.is_absolute():
            output = (paths.repo_root / output).resolve()
        return output

    return paths.profiles_dir / f"{orthogonal_path.stem}.csv"


def _iter_orthogonal_groups(orthogonal_path: Path, batch_size: int = ORTHOGONAL_READ_BATCH_SIZE):
    """Read ordered GeoPackage rows in bounded batches, retaining a split reach."""
    if batch_size < 1:
        raise ValueError("batch_size must be positive")

    with sqlite3.connect(f"file:{orthogonal_path}?mode=ro", uri=True) as connection:
        layer = connection.execute(
            "SELECT table_name FROM gpkg_contents WHERE data_type = 'features'"
        ).fetchone()[0]
        first_fid, last_fid = connection.execute(
            f'SELECT MIN(fid), MAX(fid) FROM "{layer}"'
        ).fetchone()

    if first_fid is None:
        return

    pending_id = None
    pending_group = None
    for start in range(first_fid, last_fid + 1, batch_size):
        stop = min(start + batch_size - 1, last_fid)
        batch = gpd.read_file(orthogonal_path, where=f"fid BETWEEN {start} AND {stop}")
        for combined_reach_id, group in batch.groupby("combined_reach_id", sort=False):
            if pending_group is not None and combined_reach_id == pending_id:
                pending_group = pd.concat([pending_group, group])
            else:
                if pending_group is not None:
                    yield pending_id, pending_group
                pending_id, pending_group = combined_reach_id, group

    if pending_group is not None:
        yield pending_id, pending_group


def sample_profiles_dataframe_from_orthogonals_file(
    orthogonal_path: str | Path,
    *,
    conf_factor: int | None = None,
    config_path: str | None = None,
    slope_samples: int = 400,
    dem_fill_value: int = -9999,
    dem_projection: str = "EPSG:4326",
    max_cross_distance: int = 30000,
) -> pd.DataFrame:
    paths = load_project_paths(config_path)
    paths.ensure_step2_dirs()

    orthogonal_path = Path(orthogonal_path).expanduser()
    if not orthogonal_path.is_absolute():
        orthogonal_path = (paths.repo_root / orthogonal_path).resolve()

    resolved_conf_factor = _resolve_conf_factor(orthogonal_path, conf_factor)

    if orthogonal_path.exists() is False:
        raise FileNotFoundError(f"Orthogonal intermediate not found: {orthogonal_path}")

    groups = _iter_orthogonal_groups(orthogonal_path)
    first_group = next(groups, None)
    if first_group is None:
        return empty_profile_frame()

    vrt_file = str(paths.fabdem_vrt)
    dem_vrt = gdal.Open(vrt_file)
    if dem_vrt is None:
        raise FileNotFoundError(
            f"FABDEM VRT not found at {vrt_file}. Run 'python -m pipeline.build_fabdem_index' first."
        )

    profile_rows = []
    try:
        for reach_number, (combined_reach_id, group) in enumerate(
            chain((first_group,), groups), start=1
        ):
            line_out = line_inn = left_right = centerline_wkt = np.nan
            try:
                centerline_wkt = _group_geometry(group, "centerline").wkt
                group_local = group.to_crs(group["reach_crs"].dropna().iloc[0])
                (
                    reach_crs,
                    centerline,
                    bend_lines,
                    bend_widths,
                    orthogonal_geometry,
                ) = _reconstruct_combined_reach(group_local)

                line_out = _serialize_line_array(orthogonal_geometry.line_out)
                line_inn = _serialize_line_array(orthogonal_geometry.line_inn)
                left_right = _serialize_scalar_array(orthogonal_geometry.left_right)

                sampling_df = gpd.GeoDataFrame({"geometry": [centerline]}, geometry="geometry", crs=reach_crs)
                elev_out, elev_inn, dist_out, dist_inn, line_slope, bend_height = sample_orthogonal_profiles(
                    centerline,
                    sampling_df,
                    bend_lines,
                    bend_widths,
                    orthogonal_geometry,
                    reach_crs,
                    dem_projection,
                    resolved_conf_factor,
                    dem_fill_value,
                    dem_vrt,
                    slope_samples=slope_samples,
                    maxCrossDistance=max_cross_distance,
                )

                profile_rows.append(
                    {
                        "combined_reach_id": combined_reach_id,
                        "centerlineWkt": centerline_wkt,
                        "elevOut": _serialize_profile_array(elev_out),
                        "elevInn": _serialize_profile_array(elev_inn),
                        "distOut": _serialize_profile_array(dist_out),
                        "distInn": _serialize_profile_array(dist_inn),
                        "lineOut": line_out,
                        "lineInn": line_inn,
                        "LROrthog": left_right,
                        "lineSlope": line_slope,
                        "bendHeight": _serialize_scalar_array(bend_height),
                        "sampling_code": "",
                    }
                )
            except Exception:
                profile_rows.append(
                    {
                        "combined_reach_id": combined_reach_id,
                        "centerlineWkt": centerline_wkt,
                        "elevOut": np.nan,
                        "elevInn": np.nan,
                        "distOut": np.nan,
                        "distInn": np.nan,
                        "lineOut": line_out,
                        "lineInn": line_inn,
                        "LROrthog": left_right,
                        "lineSlope": np.nan,
                        "bendHeight": np.nan,
                        "sampling_code": "4",
                    }
                )
            finally:
                # Bound cyclic garbage without scanning the whole heap twice per reach.
                if reach_number % 32 == 0:
                    gc.collect()
                if reach_number % 250 == 0:
                    rss_mb = _current_rss_mb()
                    memory_text = f", RSS {rss_mb:.0f} MB" if rss_mb is not None else ""
                    print(
                        f"Sampled profiles for {reach_number} combined reaches{memory_text}",
                        flush=True,
                    )
    finally:
        dem_vrt = None
        gc.collect()

    return pd.DataFrame(profile_rows, columns=PROFILE_COLUMNS)


def sample_profiles_from_orthogonals_file(
    orthogonal_path: str | Path,
    *,
    conf_factor: int | None = None,
    config_path: str | None = None,
    output_path: str | None = None,
    slope_samples: int = 400,
    dem_fill_value: int = -9999,
    dem_projection: str = "EPSG:4326",
    max_cross_distance: int = 30000,
) -> Path:
    paths = load_project_paths(config_path)
    paths.ensure_step2_dirs()

    orthogonal_path = Path(orthogonal_path).expanduser()
    if not orthogonal_path.is_absolute():
        orthogonal_path = (paths.repo_root / orthogonal_path).resolve()

    profile_path = _profile_output_path(paths, orthogonal_path, output_path)
    profile_path.parent.mkdir(parents=True, exist_ok=True)

    df_profiles = sample_profiles_dataframe_from_orthogonals_file(
        orthogonal_path,
        conf_factor=conf_factor,
        config_path=config_path,
        slope_samples=slope_samples,
        dem_fill_value=dem_fill_value,
        dem_projection=dem_projection,
        max_cross_distance=max_cross_distance,
    )
    df_profiles.to_csv(profile_path, index=False)
    return profile_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 2 DEM-sampling stage: sample profile values from saved orthogonal intermediates."
    )
    parser.add_argument(
        "orthogonals",
        nargs="+",
        help="One or more Step 2 orthogonal intermediate .gpkg files.",
    )
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--conf-factor",
        type=int,
        help="Cross-distance factor. If omitted, infer it from the orthogonal file name.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    for orthogonal_file in args.orthogonals:
        sample_profiles_from_orthogonals_file(
            orthogonal_file,
            conf_factor=args.conf_factor,
            config_path=args.config,
        )
