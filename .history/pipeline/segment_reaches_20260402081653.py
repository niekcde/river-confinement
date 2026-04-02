#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse
import re
from multiprocessing import Pool
from pathlib import Path

from .paths import DEFAULT_CONTINENTS, load_project_paths


SEGMENTED_OUTPUT_PATTERN = re.compile(r"^(?P<continent>[a-z]{2})_(?P<file_id>\d{2})_")

def _extract_segmented_parts(path: Path) -> tuple[str, str]:
    match = SEGMENTED_OUTPUT_PATTERN.match(path.name.lower())
    if match is None:
        raise ValueError(f"Unexpected segmented-reach filename: {path.name}")
    return match.group("continent"), match.group("file_id")


def _discover_input_pair(paths, continent: str) -> tuple[Path, Path]:
    print('debugging')
    print(paths.swot_vector_dir)
    print(paths.swot_vector_dir.glob(f"{continent}*17*gpkg"))
    vector_files = sorted(paths.swot_vector_dir.glob(f"{continent}*17*gpkg"))
    node_files = sorted(paths.swot_nodes_dir.glob(f"{continent}*17*gpkg"))

    if len(vector_files) == 0:
        raise FileNotFoundError(
            f"No Step 1 reach files found for continent '{continent}' in {paths.swot_vector_dir}"
        )
    if len(node_files) == 0:
        raise FileNotFoundError(
            f"No Step 1 node files found for continent '{continent}' in {paths.swot_nodes_dir}"
        )

    if len(vector_files) != 1 or len(node_files) != 1:
        raise ValueError(
            f"Step 1 currently expects exactly one reach file and one node file per continent. "
            f"Found {len(vector_files)} reach files and {len(node_files)} node files for '{continent}'. "
            "The downstream pipeline still assumes Step 1 outputs are named like '{continent}_{group}_...'."
        )

    return vector_files[0], node_files[0]


def _remove_existing_outputs(paths, continent: str) -> None:
    for output_file in paths.new_segments_vector_dir.glob(f"{continent}_*_reach_new_segments.gpkg"):
        output_file.unlink()
    for output_file in paths.new_segments_node_dir.glob(f"{continent}_*_node_new_segments.gpkg"):
        output_file.unlink()


def segment_continent(continent: str, config_path: str | None = None, min_reach_len_factor: int = 4 * 12, overwrite: bool = True) -> None:
    import geopandas as gpd

    from .reach_definition import new_reach_definition

    paths = load_project_paths(config_path)
    reach_file, node_file = _discover_input_pair(paths, continent)

    if overwrite:
        _remove_existing_outputs(paths, continent)

    print(f"[Step 1] segment {continent}")
    df_in = gpd.read_file(reach_file)
    df_node_in = gpd.read_file(node_file)
    new_reach_definition(df_in, df_node_in, min_reach_len_factor, paths.results_root, continent, save=True)


def build_step1_reference_tables(config_path: str | None = None) -> None:
    import geopandas as gpd
    import pandas as pd

    from .reach_definition import create_continent_new_reach
    from .support import SWORD_stats, file_sorting, smooth_factor

    paths = load_project_paths(config_path)
    files = sorted(paths.new_segments_vector_dir.glob("??_??_*.gpkg"))

    if len(files) == 0:
        raise FileNotFoundError(
            f"No segmented reach files found in {paths.new_segments_vector_dir}. "
            "Run Step 1 segmentation before building the prerequisite tables."
        )

    frames = []
    available_continents = set()
    for file_path in files:
        file_cont, file_id = _extract_segmented_parts(file_path)
        available_continents.add(file_cont)

        df = gpd.read_file(file_path)
        df["file_cont"] = file_cont
        df["file_id"] = file_id
        frames.append(df)

    df_total = pd.concat(frames, ignore_index=True)
    df_total["file"] = df_total["file_cont"] + "_" + df_total["file_id"]

    create_continent_new_reach(sorted(available_continents), paths.results_root)
    SWORD_stats(df_total, paths.results_root)
    smooth_factor(df_total, paths.results_root)
    file_sorting(df_total, paths.results_root)


def run_step1_segmentation(
    continents: list[str] | tuple[str, ...] | None = None,
    *,
    config_path: str | None = None,
    workers: int = 1,
    min_reach_len_factor: int = 4 * 12,
    overwrite: bool = True,
) -> None:
    paths = load_project_paths(config_path)
    paths.ensure_step1_dirs()
    paths.validate_step1_inputs()

    continents = list(continents or DEFAULT_CONTINENTS)
    task_args = [
        (continent, str(paths.config_path) if paths.config_path is not None else config_path, min_reach_len_factor, overwrite)
        for continent in continents
    ]

    if workers > 1 and len(task_args) > 1:
        with Pool(min(workers, len(task_args))) as pool:
            pool.starmap(segment_continent, task_args)
    else:
        for args in task_args:
            segment_continent(*args)

    build_step1_reference_tables(str(paths.config_path) if paths.config_path is not None else config_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Step 1 entrypoint: build segmented reaches and prerequisite tables."
    )
    parser.add_argument(
        "continents",
        nargs="*",
        default=list(DEFAULT_CONTINENTS),
        help="Continent codes to process. Default: af eu na sa as oc",
    )
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of continent-level worker processes to use.",
    )
    parser.add_argument(
        "--min-reach-len-factor",
        type=int,
        default=4 * 12,
        help="Minimum reach-length factor passed into new_reach_definition(...).",
    )
    parser.add_argument(
        "--keep-existing",
        action="store_true",
        help="Keep existing Step 1 outputs for a file code instead of removing and rebuilding them first.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_step1_segmentation(
        args.continents,
        config_path=args.config,
        workers=args.workers,
        min_reach_len_factor=args.min_reach_len_factor,
        overwrite=not args.keep_existing,
    )
