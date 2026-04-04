#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse

from .step2 import build_step2_orthogonals_for_file, run_geometry_cli


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 2 geometry-only entrypoint: build orthogonal intermediates without DEM sampling."
    )
    parser.add_argument("continent", nargs="?", help="Continent code to process, for example 'af' or 'oc'.")
    parser.add_argument("processors", nargs="?", type=int, help="Number of worker processes to use.")
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--vector-file",
        help="Run the geometry stage for one specific Step 1 vector file instead of a whole continent batch.",
    )
    parser.add_argument(
        "--conf-factor",
        type=int,
        default=50,
        help="Cross-distance factor used in the current Step 2 run. Default: 50.",
    )
    args = parser.parse_args()

    if args.vector_file:
        if args.continent is not None or args.processors is not None:
            parser.error("Use either 'continent processors' or '--vector-file', not both.")
    elif args.continent is None or args.processors is None:
        parser.error("Provide either 'continent processors' or '--vector-file'.")

    return args


if __name__ == "__main__":
    args = parse_args()
    if args.vector_file:
        build_step2_orthogonals_for_file(
            [args.vector_file, args.conf_factor, args.config]
        )
    else:
        run_geometry_cli(
            continent_input=args.continent,
            number_of_processors=args.processors,
            config_path=args.config,
            conf_factor=args.conf_factor,
        )
