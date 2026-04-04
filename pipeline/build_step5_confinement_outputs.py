#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

from .bend_io import read_bend_table
from .paths import load_project_paths
from .run_confinement_values import calc_confinement_values


def _format_conf_factor(conf_factor):
    if conf_factor < 10:
        return f"0{conf_factor}"
    return f"{conf_factor}"


def _format_height_factor(height_factor):
    if isinstance(height_factor, float) and height_factor.is_integer():
        height_factor = int(height_factor)
    hf_save = f"{height_factor}"
    if height_factor < 10:
        hf_save = f"0{height_factor}"
    return hf_save


def _resolve_path(paths, file_path):
    path = Path(file_path).expanduser()
    if path.is_absolute() is False:
        path = (paths.repo_root / path).resolve()
    return path


def _discover_step5_inputs(paths, continent_input=None, conf_factor=50, input_file=None):
    if input_file is not None:
        resolved = _resolve_path(paths, input_file)
        if resolved.exists() is False:
            raise FileNotFoundError(f"Step 5 input file not found: {resolved}")
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


def _resolve_step5_output_dirs(paths):
    paths.ensure_step5_dirs()
    return paths.single_values_dir, paths.reach_averaged_dir


def _step5_output_files(paths, file_stem, height_factor):
    single_values_dir, reach_averaged_dir = _resolve_step5_output_dirs(paths)
    hf_save = _format_height_factor(height_factor)
    return {
        "reach_averaged": reach_averaged_dir / f"{file_stem}_{hf_save}.gpkg",
        "bend_gpkg": single_values_dir / f"{file_stem}_{hf_save}_conf.gpkg",
        "bend_nc": single_values_dir / f"{file_stem}_{hf_save}_conf.nc",
    }


def _step5_worker(args):
    input_file, height_factor, conf_factor, config_path = args
    input_path = Path(input_file)
    if input_path.suffix.lower() == ".csv":
        df = pd.read_csv(input_path)
        open_separate = True
    else:
        df = read_bend_table(input_path)
        open_separate = False
    calc_confinement_values(
        df,
        input_path.stem,
        False,
        open_separate,
        conf_factor,
        height_factor,
        config_path=config_path,
    )
    paths = load_project_paths(config_path)
    return _step5_output_files(paths, input_path.stem, height_factor)


def run_step5_confinement_outputs(
    continent_input=None,
    *,
    config_path=None,
    workers=1,
    conf_factor=50,
    height_factor=2,
    input_file=None,
    input_files=None,
    overwrite=True,
):
    paths = load_project_paths(config_path)
    _resolve_step5_output_dirs(paths)

    if input_files is None:
        resolved_inputs = _discover_step5_inputs(
            paths,
            continent_input=continent_input,
            conf_factor=conf_factor,
            input_file=input_file,
        )
    else:
        resolved_inputs = [_resolve_path(paths, file_path) for file_path in input_files]
        if len(resolved_inputs) == 0:
            raise FileNotFoundError("No Step 5 input files were provided.")
        missing = [str(path) for path in resolved_inputs if path.exists() is False]
        if missing:
            raise FileNotFoundError(
                "Some Step 5 input files were not found:\n" + "\n".join(missing)
            )

    outputs = []
    tasks = []
    resolved_config = str(paths.config_path) if paths.config_path is not None else config_path

    for input_path in resolved_inputs:
        file_outputs = _step5_output_files(paths, input_path.stem, height_factor)
        if overwrite:
            for output_file in file_outputs.values():
                if output_file.exists():
                    output_file.unlink()
        elif all(output_file.exists() for output_file in file_outputs.values()):
            outputs.append(file_outputs)
            continue

        tasks.append((str(input_path), height_factor, conf_factor, resolved_config))

    if workers == 1 or len(tasks) <= 1:
        for task in tasks:
            outputs.append(_step5_worker(task))
    else:
        with Pool(min(workers, len(tasks))) as pool:
            for output in pool.imap(_step5_worker, tasks):
                outputs.append(output)

    return outputs


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Step 5 entrypoint: compute confinement outputs from bend-level Step 2 tables."
    )
    parser.add_argument("continent", nargs="?", help="Continent code to process, for example 'af' or 'oc'.")
    parser.add_argument("processors", nargs="?", type=int, help="Number of worker processes to use.")
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--conf-factor",
        type=int,
        default=50,
        help="Cross-distance factor used in the current Step 5 run. Default: 50.",
    )
    parser.add_argument(
        "--height-factor",
        type=float,
        default=2,
        help="Height factor for the current Step 5 run. Default: 2.",
    )
    parser.add_argument(
        "--input-file",
        help="Run Step 5 for one specific bend-level Step 2 file instead of a whole continent batch.",
    )
    parser.add_argument(
        "--keep-existing",
        action="store_true",
        help="Keep existing Step 5 outputs instead of rebuilding them first.",
    )
    args = parser.parse_args(argv)

    if args.input_file:
        if args.continent is not None or args.processors is not None:
            parser.error("Use either 'continent processors' or '--input-file', not both.")
    elif args.continent is None:
        parser.error("Provide either 'continent processors' or '--input-file'.")
    elif args.processors is None:
        parser.error("Missing required positional argument: processors")

    return args


def main_cli(argv=None):
    args = parse_args(argv)
    workers = 1 if args.input_file else args.processors
    run_step5_confinement_outputs(
        continent_input=args.continent,
        config_path=args.config,
        workers=workers,
        conf_factor=args.conf_factor,
        height_factor=args.height_factor,
        input_file=args.input_file,
        overwrite=not args.keep_existing,
    )


if __name__ == "__main__":
    main_cli()
