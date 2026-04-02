#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse

from .run_confinement_values import concat_nc_conf_files, concat_reachAveraged


def run_step6_aggregates(*, cross_factor=50, height_factor=2, config_path=None, directory=None):
    global_nc = concat_nc_conf_files(
        directory=directory,
        cross=cross_factor,
        hf=height_factor,
        config_path=config_path,
    )
    reach_outputs = concat_reachAveraged(
        directory=directory,
        cross=cross_factor,
        hf=height_factor,
        config_path=config_path,
    )
    return {
        "global_nc": global_nc,
        "reach_averaged": reach_outputs,
    }


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Step 6 entrypoint: aggregate Step 5 outputs for one height factor."
    )
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--cross-factor",
        type=int,
        default=50,
        help="Cross-distance factor used in the current Step 6 run. Default: 50.",
    )
    parser.add_argument(
        "--height-factor",
        type=float,
        required=True,
        help="Height factor to aggregate, for example 2, 1.5, or 0.5.",
    )
    return parser.parse_args(argv)


def main_cli(argv=None):
    args = parse_args(argv)
    run_step6_aggregates(
        cross_factor=args.cross_factor,
        height_factor=args.height_factor,
        config_path=args.config,
    )


if __name__ == "__main__":
    main_cli()
