#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse
from pathlib import Path

import pandas as pd

from .bend_io import read_bend_table
from .paths import load_project_paths
from .support import confinement_factor_single_values


def _format_conf_factor(conf_factor):
    if conf_factor < 10:
        return f"0{conf_factor}"
    return f"{conf_factor}"


def _resolve_path(paths, file_path):
    path = Path(file_path).expanduser()
    if path.is_absolute() is False:
        path = (paths.repo_root / path).resolve()
    return path


def _discover_step4_inputs(paths, continent_input=None, conf_factor=50, input_file=None):
    if input_file is not None:
        resolved = _resolve_path(paths, input_file)
        if resolved.exists() is False:
            raise FileNotFoundError(f"Step 4 input file not found: {resolved}")
        return [resolved]

    conf_factor_str = _format_conf_factor(conf_factor)
    if continent_input is None:
        files = sorted(paths.bends_dir.glob(f"??_??_{conf_factor_str}.parquet"))
    else:
        files = sorted(paths.bends_dir.glob(f"{continent_input}_??_{conf_factor_str}.parquet"))

    if len(files) == 0:
        if continent_input is None:
            raise FileNotFoundError(
                f"No bend-level Step 2 files found in {paths.bends_dir} for conf_factor {conf_factor_str}."
            )
        raise FileNotFoundError(
            f"No bend-level Step 2 files found for continent '{continent_input}' in {paths.bends_dir} "
            f"for conf_factor {conf_factor_str}."
        )

    return files


def _resolve_output_file(paths, conf_factor=50, output_path=None):
    if output_path is None:
        return paths.confinement_factor_file_for(conf_factor)

    output_file = Path(output_path).expanduser()
    if output_file.is_absolute() is False:
        output_file = (paths.repo_root / output_file).resolve()
    return output_file


def build_step4_confinement_factor(
    continent_input=None,
    *,
    config_path=None,
    conf_factor=50,
    input_file=None,
    input_files=None,
    output_path=None,
    overwrite=True,
    width_column="bendWidths",
    y1=50,
    y2=10,
):
    paths = load_project_paths(config_path)
    paths.ensure_step4_dirs()

    if input_files is None:
        resolved_inputs = _discover_step4_inputs(
            paths,
            continent_input=continent_input,
            conf_factor=conf_factor,
            input_file=input_file,
        )
    else:
        resolved_inputs = [_resolve_path(paths, file_path) for file_path in input_files]
        if len(resolved_inputs) == 0:
            raise FileNotFoundError("No bend-level Step 2 input files were provided for Step 4.")
        missing = [str(path) for path in resolved_inputs if path.exists() is False]
        if missing:
            raise FileNotFoundError(
                "Some Step 4 input files were not found:\n" + "\n".join(missing)
            )

    output_file = _resolve_output_file(paths, conf_factor, output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.exists():
        if overwrite:
            output_file.unlink()
        else:
            return output_file

    width_frames = [read_bend_table(input_path)[[width_column]] for input_path in resolved_inputs]
    if len(width_frames) == 0:
        raise FileNotFoundError("No width values were loaded for Step 4.")

    df_widths = pd.concat(width_frames, ignore_index=True)
    if df_widths.empty:
        raise ValueError("Step 4 loaded zero bend-width rows from the available bend-level Step 2 files.")

    df_factor = confinement_factor_single_values(df_widths, width_column, y1, y2)
    df_factor.to_csv(output_file, index=False)
    return output_file


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Step 4 entrypoint: build the confinement-factor lookup table from bend-level Step 2 outputs."
    )
    parser.add_argument("continent", nargs="?", help="Continent code to process, for example 'af' or 'oc'.")
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--conf-factor",
        type=int,
        default=50,
        help="Cross-distance factor used in the current Step 4 run. Default: 50.",
    )
    parser.add_argument(
        "--input-file",
        help="Run Step 4 for one specific bend-level Step 2 file instead of a whole continent batch.",
    )
    parser.add_argument(
        "--output-file",
        help="Optional custom output path for the Step 4 confinement-factor table.",
    )
    parser.add_argument(
        "--keep-existing",
        action="store_true",
        help="Keep an existing output file instead of rebuilding it first.",
    )
    args = parser.parse_args(argv)
    if args.input_file and args.continent is not None:
        parser.error("Use either 'continent' or '--input-file', not both.")
    return args


def main_cli(argv=None):
    args = parse_args(argv)
    build_step4_confinement_factor(
        continent_input=args.continent,
        config_path=args.config,
        conf_factor=args.conf_factor,
        input_file=args.input_file,
        output_path=args.output_file,
        overwrite=not args.keep_existing,
    )


if __name__ == "__main__":
    main_cli()
