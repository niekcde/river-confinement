#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse
import gc
from multiprocessing import Pool
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from .dem import find_dem_FAB, find_dem_bounds_FAB
from .paths import load_project_paths
from .support import adjust_new_segments, check_memory


def divide_dataframe_in_equal_parts(df, max_size=1000):
    new_rows = []

    for _, row in df.iterrows():
        value = row["size"]
        parts = int(np.ceil(value / max_size))
        step = value // parts
        remainder = value % parts

        start = 0
        for i in range(parts):
            end = start + step + (remainder if i == (parts - 1) else 0)
            new_rows.append(
                {"file": row["file"], "values": value, "start": start, "end": end, "order": i}
            )
            start = end

    new_df = pd.DataFrame(new_rows)
    return df.merge(new_df, on="file", how="left")


def _select_raster_dirs(paths):
    raster_dir = paths.input_created_root / "FAB_dem_reach"
    dataframe_dir = paths.input_created_root / "FAB_dem_reach_dataframe"
    raster_dir.mkdir(parents=True, exist_ok=True)
    dataframe_dir.mkdir(parents=True, exist_ok=True)
    return raster_dir, dataframe_dir


def save_raster_reach(dfReach, reachID, reachCRS, cont, bufferSize, dfDemBounds, *, raster_dir):
    raster, dfDemRows = find_dem_FAB(dfReach, bufferSize, dfDemBounds, reachCRS, "EPSG:4326")
    if dfDemRows.shape[0] == 0:
        print(dfReach.reach_id.values)
        rasterPath = ""
        rasterType = "None"
        rasterBounds = ""
        rasterCRS = ""
    elif dfDemRows.shape[0] == 1:
        rasterPath = dfDemRows.iloc[0]["id"]
        rasterType = "single"
        rasterBounds = str(raster.rio.bounds())
        rasterCRS = str(raster.rio.crs)
    else:
        raster_path = raster_dir / f"{cont}_{int(reachID)}_DEM.tif"
        raster.rio.to_raster(raster_path)
        rasterPath = str(raster_path)
        rasterType = "multi"
        rasterBounds = str(raster.rio.bounds())
        rasterCRS = str(raster.rio.crs)

    dfDEM = pd.DataFrame(
        {
            "combined_reach_id": reachID,
            "Continent": cont,
            "bounds": rasterBounds,
            "bounds_crs": rasterCRS,
            "rasterPath": rasterPath,
            "rasterType": rasterType,
        },
        index=[0],
    )
    raster.close()
    del raster
    gc.collect()
    return dfDEM


def run_save_raster_reach(df, cont, contOrder, bufferSize, dfDemBounds, *, paths):
    raster_dir, dataframe_dir = _select_raster_dirs(paths)
    dfDEM = 0
    ids = df.loc[df["include_flag"] == "0", "combined_reach_id"].unique()

    for rid in ids:
        dfReach = df[df["combined_reach_id"] == rid].copy()
        groupedCRS = dfReach.groupby("localCRS", as_index=False).size()
        reachCRS = groupedCRS.loc[
            groupedCRS["size"] == groupedCRS["size"].max(),
            "localCRS",
        ].iloc[0]
        dfReach = dfReach.to_crs(reachCRS)

        dfDEMTemp = save_raster_reach(
            dfReach,
            rid,
            reachCRS,
            cont,
            dfReach["combined_reach_max_width"].iloc[0] * bufferSize,
            dfDemBounds,
            raster_dir=raster_dir,
        )
        if isinstance(dfDEM, int):
            dfDEM = dfDEMTemp
        else:
            dfDEM = pd.concat([dfDEM, dfDEMTemp])

    output_file = dataframe_dir / f"{cont}_{contOrder}_{bufferSize}.csv"
    dfDEM.to_csv(output_file, index=False)
    return output_file


def _continent_from_vector_path(file_path):
    return Path(file_path).name[:2]


def multi_save_raster_reach(multiInput):
    filePath, start, end, order, bufferSize, config_path = multiInput
    paths = load_project_paths(config_path)
    dfDemBounds = find_dem_bounds_FAB(
        demCRS="EPSG:4326",
        fabdem_dir=paths.fabdem_dir,
        dem_boundary_file=paths.fabdem_bounds,
    )

    df = gpd.read_file(filePath)
    df = adjust_new_segments(df)
    size = df.loc[df["include_flag"] == "0"].shape[0]
    if end == size:
        dfFilter = df.iloc[start::]
    else:
        dfFilter = df.iloc[start:end]

    cont = _continent_from_vector_path(filePath)
    print(cont, order, check_memory())
    output_file = run_save_raster_reach(
        dfFilter,
        cont,
        order,
        bufferSize,
        dfDemBounds,
        paths=paths,
    )
    print(cont, order, check_memory())
    print()
    del df, dfFilter
    return output_file


def run_select_raster_chunks(*, config_path=None, workers=1, buffer_size=50, limit_tasks=None, max_size=1000):
    paths = load_project_paths(config_path)
    _, dataframe_dir = _select_raster_dirs(paths)

    sortFiles = pd.read_csv(paths.reference_tables_dir / "file_sorting.csv").sort_values("size")
    newSortFiles = divide_dataframe_in_equal_parts(sortFiles, max_size=max_size)
    multiInput = newSortFiles[["filePath", "start", "end", "order"]].values.tolist()
    if limit_tasks is not None:
        multiInput = multiInput[:limit_tasks]

    resolved_config = str(paths.config_path) if paths.config_path is not None else config_path
    tasks = [tuple(list(row) + [buffer_size, resolved_config]) for row in multiInput]

    for csv_file in dataframe_dir.glob("*"):
        csv_file.unlink()

    print("start Multiprocess", len(tasks))
    print(tasks)
    if workers == 1 or len(tasks) <= 1:
        outputs = [multi_save_raster_reach(task) for task in tasks]
    else:
        with Pool(min(workers, len(tasks))) as p:
            outputs = list(p.imap(multi_save_raster_reach, tasks))

    print("End Multiprocess")
    files = sorted(dataframe_dir.glob("*.csv"))
    if len(files) == 0:
        raise FileNotFoundError(
            f"No raster-reach dataframe chunks were produced in {dataframe_dir}."
        )

    frames = [pd.read_csv(f) for f in files]
    dfR = pd.concat(frames, ignore_index=True)
    final_output = dataframe_dir / "FAB_dem_reach.csv"
    dfR.to_csv(final_output, index=False)
    print("code Finished")
    return {
        "chunk_outputs": outputs,
        "final_output": final_output,
    }


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Utility entrypoint: build FABDEM reach-raster chunk tables from segmented Step 1 outputs."
    )
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of worker processes to use. Default: 1.",
    )
    parser.add_argument(
        "--buffer-size",
        type=float,
        default=50,
        help="Buffer multiplier applied to combined reach max width. Default: 50.",
    )
    parser.add_argument(
        "--limit-tasks",
        type=int,
        help="Optional limit on the number of file chunks to process, useful for debugging.",
    )
    parser.add_argument(
        "--max-size",
        type=int,
        default=1000,
        help="Maximum target chunk size used when splitting file_sorting.csv groups. Default: 1000.",
    )
    return parser.parse_args(argv)


def main_cli(argv=None):
    args = parse_args(argv)
    run_select_raster_chunks(
        config_path=args.config,
        workers=args.workers,
        buffer_size=args.buffer_size,
        limit_tasks=args.limit_tasks,
        max_size=args.max_size,
    )


if __name__ == "__main__":
    main_cli()
