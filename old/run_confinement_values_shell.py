#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse
from datetime import datetime as dt

from .build_step3_single_values import run_step3_single_values
from .build_step4_confinement_factor import build_step4_confinement_factor
from .build_step5_confinement_outputs import run_step5_confinement_outputs
from .build_step6_aggregates import run_step6_aggregates
from .paths import load_project_paths


DEFAULT_HEIGHT_FACTORS = [2, 0.5, 1, 1.5, 3, 4, 6, 8, 10, 15]
DEFAULT_CON_FACTOR = [50, 10]


def run_legacy_confinement_workflow(
    continent,
    *,
    workers=10,
    config_path=None,
    cross_factor=50,
    height_factors=None,
    con_factor=None,
    create_new_single_values=True,
    create_new_factor=True,
):
    paths = load_project_paths(config_path)
    height_factors = DEFAULT_HEIGHT_FACTORS if height_factors is None else list(height_factors)
    con_factor = DEFAULT_CON_FACTOR if con_factor is None else list(con_factor)

    dt1 = dt.now()
    print(f"start Code: {dt1}")

    all_result_files = None
    if create_new_single_values:
        print("Start transform results in single row values")
        step3_outputs = run_step3_single_values(
            continent_input=continent,
            workers=workers,
            config_path=config_path,
            conf_factor=cross_factor,
        )
        all_result_files = [str(path) for path in step3_outputs]
        dt2 = dt.now()
        print(f"Open_to_single_apex Finished: {dt2 - dt1}")
    else:
        conf_token = f"{cross_factor:02d}"
        all_result_files = [
            str(path)
            for path in sorted(paths.single_values_dir.glob(f"{continent}_??_{conf_token}.csv"))
        ]
        if len(all_result_files) == 0:
            raise FileNotFoundError(
                f"No Step 3 files found for continent '{continent}' in {paths.single_values_dir} "
                f"for conf_factor {conf_token}."
            )
        dt2 = dt1

    if create_new_factor:
        print("Start calc confinement factor")
        build_step4_confinement_factor(
            input_files=all_result_files,
            config_path=config_path,
            conf_factor=cross_factor,
            y1=con_factor[0],
            y2=con_factor[1],
        )
        dt3 = dt.now()
        print(f"Confinement factor Finished: {dt3 - dt2}")
    else:
        dt3 = dt2

    print("Start calc confinement values")
    for hf in height_factors:
        run_step5_confinement_outputs(
            input_files=all_result_files,
            workers=workers,
            conf_factor=cross_factor,
            height_factor=hf,
            config_path=config_path,
        )
        run_step6_aggregates(
            cross_factor=cross_factor,
            height_factor=hf,
            config_path=config_path,
        )
        print(f"run confinement heightfactor {hf} finished")

    print(f"run confinement finished Finished: {dt.now() - dt1}")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Legacy wrapper for Steps 3 to 6 using the canonical stage entrypoints."
    )
    parser.add_argument("continent", help="Continent code to process, for example 'af' or 'oc'.")
    parser.add_argument(
        "--workers",
        type=int,
        default=10,
        help="Number of worker processes to use for Steps 3 and 5. Default: 10.",
    )
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--cross-factor",
        type=int,
        default=50,
        help="Cross-distance factor used in the current run. Default: 50.",
    )
    parser.add_argument(
        "--height-factors",
        nargs="*",
        type=float,
        default=list(DEFAULT_HEIGHT_FACTORS),
        help="Height factors to process. Default: 2 0.5 1 1.5 3 4 6 8 10 15.",
    )
    parser.add_argument(
        "--con-factor",
        nargs=2,
        type=float,
        default=list(DEFAULT_CON_FACTOR),
        metavar=("Y1", "Y2"),
        help="Confinement-factor curve control points. Default: 50 10.",
    )
    parser.add_argument(
        "--skip-step3",
        action="store_true",
        help="Reuse existing Step 3 outputs instead of rebuilding them.",
    )
    parser.add_argument(
        "--skip-step4",
        action="store_true",
        help="Reuse the existing Step 4 confinement-factor table instead of rebuilding it.",
    )
    return parser.parse_args(argv)


def main_cli(argv=None):
    args = parse_args(argv)
    run_legacy_confinement_workflow(
        args.continent,
        workers=args.workers,
        config_path=args.config,
        cross_factor=args.cross_factor,
        height_factors=args.height_factors,
        con_factor=args.con_factor,
        create_new_single_values=not args.skip_step3,
        create_new_factor=not args.skip_step4,
    )


if __name__ == "__main__":
    main_cli()
