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

from .paths import load_project_paths
from .run_confinement_values import create_apex_val_dataframe


def _format_conf_factor(conf_factor):
    if conf_factor < 10:
        return f"0{conf_factor}"
    return f"{conf_factor}"


def _step3_worker(args):
    input_file, output_file, config_path = args
    create_apex_val_dataframe(
        input_file,
        output_path=output_file,
        config_path=config_path,
    )
    return output_file


def _resolve_input_file(paths, input_file):
    path = Path(input_file).expanduser()
    if path.is_absolute() is False:
        path = (paths.repo_root / path).resolve()
    return path


def _discover_step3_inputs(paths, continent_input=None, conf_factor=50, input_file=None):
    if input_file is not None:
        resolved = _resolve_input_file(paths, input_file)
        if resolved.exists() is False:
            raise FileNotFoundError(f"Step 3 input file not found: {resolved}")
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


def run_step3_single_values(
    continent_input=None,
    *,
    config_path=None,
    workers=1,
    conf_factor=50,
    input_file=None,
    overwrite=True,
):
    paths = load_project_paths(config_path)
    paths.ensure_step3_dirs()

    input_files = _discover_step3_inputs(
        paths,
        continent_input=continent_input,
        conf_factor=conf_factor,
        input_file=input_file,
    )
    output_files = [paths.single_values_dir / f"{input_path.stem}.csv" for input_path in input_files]

    if overwrite:
        for output_file in output_files:
            if output_file.exists():
                output_file.unlink()

    resolved_config = str(paths.config_path) if paths.config_path is not None else config_path
    tasks = [(str(input_file), str(output_file), resolved_config) for input_file, output_file in zip(input_files, output_files)]

    if workers == 1 or len(tasks) == 1:
        for task in tasks:
            _step3_worker(task)
    else:
        with Pool(min(workers, len(tasks))) as pool:
            for _ in pool.imap(_step3_worker, tasks):
                pass

    return output_files


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Step 3 compatibility entrypoint: export bend-level Step 2 tables to the legacy single_values CSV format."
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
        help="Cross-distance factor used in the current Step 3 run. Default: 50.",
    )
    parser.add_argument(
        "--input-file",
        help="Run Step 3 for one specific bend-level Step 2 file instead of a whole continent batch.",
    )
    parser.add_argument(
        "--keep-existing",
        action="store_true",
        help="Keep existing Step 3 outputs instead of rebuilding them first.",
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
    if args.input_file:
        workers = 1
    else:
        workers = args.processors

    run_step3_single_values(
        continent_input=args.continent,
        config_path=args.config,
        workers=workers,
        conf_factor=args.conf_factor,
        input_file=args.input_file,
        overwrite=not args.keep_existing,
    )


if __name__ == "__main__":
    main_cli()
