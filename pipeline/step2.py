#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# """
# Created on Tue Apr 30 10:36:11 2024
# @author: niekcollotdescury
# """

import glob
import os
import warnings
from datetime import datetime as dt
from multiprocessing import Pool

warnings.filterwarnings("ignore")

from .paths import load_project_paths

os.environ["MALLOC_MMAP_MAX_"] = "40960"


def ensure_main_dirs(paths):
    for path in (
        paths.input_created_root,
        paths.input_created_dem_dir,
        paths.results_root,
        paths.reference_tables_dir,
        paths.new_segments_root,
        paths.new_segments_node_dir,
        paths.new_segments_vector_dir,
        paths.orthogonals_dir,
        paths.profiles_dir,
        paths.all_dir,
        paths.centerline_dir,
        paths.cycles_dir,
    ):
        path.mkdir(parents=True, exist_ok=True)


def _step2_file_parts(file, conf_factor):
    cont_name = file[-29:-27]
    file_number = file[-26:-24]
    if conf_factor < 10:
        csd = f"0{conf_factor}"
    else:
        csd = f"{conf_factor}"
    file_stem = f"{cont_name}_{file_number}_{csd}"
    return cont_name, file_number, file_stem


def _build_step2_geometry_data(file, conf_factor, config_path):
    import geopandas as gpd
    import numpy as np
    import pandas as pd
    from scipy.stats import trim_mean

    from .connect_geometries import merge_centerlines
    from .get_orthogonals import build_orthogonal_lines, build_orthogonal_stage_frame
    from .inflection_points import inflection_points_curve
    from .smoothing import SG_smoothing
    from .support import adjust_new_segments, node_position

    paths = load_project_paths(config_path)
    ensure_main_dirs(paths)

    cont_name, file_number, file_stem = _step2_file_parts(file, conf_factor)
    print(f"Start: {cont_name}_{file_number}. Time: {dt.now()}")
    code_failure = 0
    dem_projection = "EPSG:4326"

    vector_file = glob.glob(str(paths.new_segments_vector_dir / f"{cont_name}_{file_number}_*.gpkg"))[0]
    node_file = glob.glob(str(paths.new_segments_node_dir / f"{cont_name}_{file_number}_*.gpkg"))[0]
    df = gpd.read_file(vector_file)
    df_node = gpd.read_file(node_file)
    print("read file finished")
    df = adjust_new_segments(df)

    included_mask = df["include_flag"] == "0"
    included_df = df.loc[included_mask].copy()
    ids = included_df["combined_reach_id"].unique()

    reach_groups = {
        combined_reach_id: df.loc[group.index].copy()
        for combined_reach_id, group in included_df.groupby("combined_reach_id", sort=False)
    }
    reach_crs_lookup = (
        included_df.groupby("combined_reach_id")["localCRS"]
        .agg(lambda values: values.value_counts().idxmax())
        .to_dict()
    )

    reach_to_combined = df[["reach_id", "combined_reach_id"]].drop_duplicates()
    df_node_with_combined = df_node.merge(
        reach_to_combined.rename(columns={"combined_reach_id": "combined_reach_id_from_reach"}),
        on="reach_id",
        how="left",
    )
    if "combined_reach_id" in df_node_with_combined.columns:
        df_node_with_combined["combined_reach_id"] = df_node_with_combined["combined_reach_id"].fillna(
            df_node_with_combined["combined_reach_id_from_reach"]
        )
        df_node_with_combined = df_node_with_combined.drop(columns=["combined_reach_id_from_reach"])
    else:
        df_node_with_combined = df_node_with_combined.rename(
            columns={"combined_reach_id_from_reach": "combined_reach_id"}
        )
    node_groups = {
        combined_reach_id: group.drop(columns=["combined_reach_id"]).copy()
        for combined_reach_id, group in df_node_with_combined.groupby("combined_reach_id", sort=False)
    }
    df_sf = pd.read_csv(paths.reference_tables_dir / "smoothingFactor.csv")
    orthogonal_frames = []

    print(f"Start for loop {len(ids)}")
    for current_id in ids:
        calculated = "0"
        df_reach = reach_groups[current_id].copy()
        df_reach_nodes = node_groups.get(current_id)
        if df_reach_nodes is None:
            print(current_id, "missing_nodes")
            continue
        df_reach_nodes = df_reach_nodes.copy()
        reach_crs = reach_crs_lookup[current_id]

        try:
            df_reach = df_reach.to_crs(reach_crs)
            df_reach_nodes = df_reach_nodes.to_crs(reach_crs)
            df.loc[df_reach.index, "combined_reach_crs"] = reach_crs

            combined_line, _, _ = merge_centerlines(df_reach, df, reach_crs)

            width_var = "width"
            factor_width = df_reach["combined_reach_width"].iloc[0]

            if df_reach_nodes[width_var].std() > (factor_width / 2):
                factor_width = trim_mean(df_reach_nodes[width_var], 0.05)
            if np.isnan(factor_width):
                factor_width = df_reach_nodes[width_var].mean()

            width_ratio = df_reach_nodes["width"] / df_reach_nodes["max_width"]
            df.loc[df_reach.index, "widthRatio"] = np.nanmean(width_ratio)

            factor_row = abs(df_sf["combined_reach_width"] - factor_width).argsort()
            smoothing_factor = df_sf.loc[factor_row, "smoothFactor"].iloc[0]
            smoothing_window = smoothing_factor * int(factor_width)

            combined_line = SG_smoothing(combined_line, smoothing_window, factor_width, id=current_id)
            if len(combined_line.coords) < 3:
                combined_line = combined_line.segmentize(1)

            df_reach_nodes = node_position(combined_line, df_reach_nodes)

        except Exception:
            calculated += "1"
            print(current_id, calculated)

        try:
            (
                sin,
                bend_sin,
                inf_p,
                inf_p_total,
                apex,
                apex_p,
                apex_po,
                ang,
                bend_lines,
                inf_lines,
                bend_widths,
                bend_max_widths,
                bend_dist_out,
                bend_len,
                bend_facc,
            ) = inflection_points_curve(combined_line, df_reach, df_reach_nodes)

            bend_sin = np.array2string(bend_sin, separator=", ")
            ang = np.array2string(ang, separator=", ")
            bend_dist_out = np.array2string(bend_dist_out, separator=", ")
            bend_facc = np.array2string(bend_facc, separator=", ")
            bend_len = np.array2string(bend_len, separator=", ")

            apex_str = np.array2string(apex, separator=", ")
            bend_widths_str = np.array2string(bend_widths, separator=", ")
            bend_max_widths_str = np.array2string(bend_max_widths, separator=", ")

            apex_p_wkt = str([geometry.wkt for geometry in apex_p])
            bend_lines_wkt = str([geometry.wkt for geometry in bend_lines])
            inf_lines_wkt = str([geometry.wkt for geometry in inf_lines])

        except Exception:
            sin, bend_sin, apex, apex_p, apex_po, ang = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
            inf_lines, bend_lines = np.nan, np.nan
            bend_widths, bend_max_widths = np.nan, np.nan
            bend_dist_out, bend_facc, bend_len = np.nan, np.nan, np.nan

            apex_str, bend_widths_str, bend_max_widths_str = np.nan, np.nan, np.nan
            apex_p_wkt, inf_lines_wkt, bend_lines_wkt = np.nan, np.nan, np.nan

            code_failure += 1
            calculated += "2"

        df.loc[df_reach.index, "sin"] = sin
        df.loc[df_reach.index, "bendSin"] = bend_sin
        df.loc[df_reach.index, "apex"] = apex_str
        df.loc[df_reach.index, "apexP"] = apex_p_wkt
        df.loc[df_reach.index, "ang"] = ang

        df.loc[df_reach.index, "infP"] = inf_lines_wkt
        df.loc[df_reach.index, "bendLines"] = bend_lines_wkt
        df.loc[df_reach.index, "bendWidths"] = bend_widths_str
        df.loc[df_reach.index, "bendMaxWidths"] = bend_max_widths_str
        df.loc[df_reach.index, "bendDistOut"] = bend_dist_out
        df.loc[df_reach.index, "bendLen"] = bend_len
        df.loc[df_reach.index, "bendFacc"] = bend_facc

        try:
            if isinstance(apex_p, np.ndarray) is True:
                orthogonal_geometry = build_orthogonal_lines(
                    combined_line,
                    apex_p,
                    apex_po,
                    apex,
                    inf_lines,
                    bend_widths,
                    conf_factor,
                )

                orthogonal_frame = build_orthogonal_stage_frame(
                    current_id,
                    reach_crs,
                    combined_line,
                    bend_lines,
                    bend_widths,
                    orthogonal_geometry,
                )
                if orthogonal_frame is None:
                    calculated += "3"
                else:
                    orthogonal_frames.append(orthogonal_frame.to_crs(dem_projection))
            else:
                calculated += "3"

        except Exception:
            code_failure += 1
            calculated += "3"

        df.loc[df_reach.index, "calculated"] = calculated

        del df_reach, df_reach_nodes
        del apex_p, inf_lines, bend_lines

    orthogonal_path = paths.orthogonals_dir / f"{file_stem}.gpkg"
    profile_path = paths.profiles_dir / f"{file_stem}.csv"

    if orthogonal_path.exists():
        orthogonal_path.unlink()

    orthogonal_written = False
    if orthogonal_frames:
        orthogonal_gdf = gpd.GeoDataFrame(
            pd.concat(orthogonal_frames, ignore_index=True),
            geometry="geometry",
            crs=dem_projection,
        )
        orthogonal_gdf.to_file(orthogonal_path, driver="GPKG")
        orthogonal_written = True

    return {
        "paths": paths,
        "cont_name": cont_name,
        "file_number": file_number,
        "file_stem": file_stem,
        "conf_factor": conf_factor,
        "config_path": config_path,
        "code_failure": code_failure,
        "df": df,
        "included_mask": included_mask,
        "orthogonal_path": orthogonal_path,
        "profile_path": profile_path,
        "orthogonal_written": orthogonal_written,
    }


def build_step2_orthogonals_for_file(multi_input):
    file, conf_factor, config_path = multi_input
    geometry_data = _build_step2_geometry_data(file, conf_factor, config_path)
    print(
        f"Finish geometry: {geometry_data['cont_name']}_{geometry_data['file_number']} "
        f"with cross slope {conf_factor}, code failures: {geometry_data['code_failure']}. Time: {dt.now()}"
    )
    return {
        "file_stem": geometry_data["file_stem"],
        "orthogonal_path": str(geometry_data["orthogonal_path"]) if geometry_data["orthogonal_written"] else None,
        "code_failure": geometry_data["code_failure"],
    }


def process_file(multi_input):
    import numpy as np
    import pandas as pd

    from .sample_step2_profiles import PROFILE_COLUMNS, empty_profile_frame, sample_profiles_from_orthogonals_file

    file, conf_factor, config_path = multi_input
    geometry_data = _build_step2_geometry_data(file, conf_factor, config_path)

    paths = geometry_data["paths"]
    df = geometry_data["df"]
    included_mask = geometry_data["included_mask"]
    profile_path = geometry_data["profile_path"]

    if profile_path.exists():
        profile_path.unlink()

    if geometry_data["orthogonal_written"]:
        sample_profiles_from_orthogonals_file(
            geometry_data["orthogonal_path"],
            conf_factor=conf_factor,
            config_path=config_path,
            output_path=str(profile_path),
        )
    else:
        empty_profile_frame().to_csv(profile_path, index=False)

    df_profiles = pd.read_csv(profile_path)
    for column in [column for column in PROFILE_COLUMNS if column not in {"combined_reach_id", "sampling_code"}]:
        if column not in df.columns:
            df[column] = np.nan

    if df_profiles.empty is False:
        df_profiles = df_profiles.drop_duplicates(subset=["combined_reach_id"]).set_index("combined_reach_id")
        included_ids = df.loc[included_mask, "combined_reach_id"]

        for column in [column for column in PROFILE_COLUMNS if column not in {"combined_reach_id", "sampling_code"}]:
            df.loc[included_mask, column] = included_ids.map(df_profiles[column])

        sampling_codes = included_ids.map(df_profiles["sampling_code"]).fillna("")
        has_sampling_code = sampling_codes != ""
        sampling_indices = sampling_codes.index[has_sampling_code]
        df.loc[sampling_indices, "calculated"] = (
            df.loc[sampling_indices, "calculated"].fillna("0")
            + sampling_codes.loc[sampling_indices]
        )

    file_name = f"{geometry_data['file_stem']}.csv"
    df.to_csv(paths.all_dir / file_name, index=False)

    print(
        f"Finish: {geometry_data['cont_name']}_{geometry_data['file_number']} with cross slope {conf_factor}, "
        f"code failures: {geometry_data['code_failure']}. Time: {dt.now()}"
    )


def step2_multi_input(continent_input, config_path=None, conf_factor=50, number_of_processors=None):
    import pandas as pd

    paths = load_project_paths(config_path)
    ensure_main_dirs(paths)

    df_files = pd.read_csv(paths.reference_tables_dir / "file_sorting.csv", index_col=0)
    df_files = df_files[df_files["file"].str.startswith(continent_input)].sort_values("size", ascending=False)
    files = df_files["filePath"].values
    if number_of_processors is None:
        print(files, continent_input)
    else:
        print(files, continent_input, number_of_processors)
    print()

    resolved_config = str(paths.config_path) if paths.config_path is not None else config_path
    return [[file, conf_factor, resolved_config] for file in files]


def run_geometry_cli(continent_input, number_of_processors, config_path=None, conf_factor=50):
    paths = load_project_paths(config_path)
    ensure_main_dirs(paths)

    remove_files = glob.glob(str(paths.orthogonals_dir / f"{continent_input}*_{conf_factor}.gpkg"))
    for remove_file in remove_files:
        os.remove(remove_file)

    multi_input = step2_multi_input(
        continent_input=continent_input,
        config_path=config_path,
        conf_factor=conf_factor,
        number_of_processors=number_of_processors,
    )

    if number_of_processors == 1:
        for item in multi_input:
            build_step2_orthogonals_for_file(item)
    else:
        with Pool(number_of_processors) as pool:
            for _ in pool.imap(build_step2_orthogonals_for_file, multi_input):
                pass

    print("geometry stage finished")


def run_results_cli(continent_input, number_of_processors, config_path=None, conf_factor=50):
    paths = load_project_paths(config_path)
    ensure_main_dirs(paths)

    remove_files = glob.glob(str(paths.all_dir / f"{continent_input}*_{conf_factor}.csv"))
    for remove_file in remove_files:
        os.remove(remove_file)

    multi_input = step2_multi_input(
        continent_input=continent_input,
        config_path=config_path,
        conf_factor=conf_factor,
        number_of_processors=number_of_processors,
    )

    if number_of_processors == 1:
        for item in multi_input:
            process_file(item)
    else:
        with Pool(number_of_processors) as pool:
            for _ in pool.imap(process_file, multi_input):
                pass

    print("code Finished")
