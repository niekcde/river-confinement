#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse

from .step2 import process_file, run_results_cli


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Step 2 entrypoint: build orthogonals, sample DEM profiles, and write bend-level outputs."
    )
    parser.add_argument("continent", nargs="?", help="Continent code to process, for example 'af' or 'eu'.")
    parser.add_argument("processors", nargs="?", type=int, help="Number of worker processes to use.")
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--vector-file",
        help="Run the full Step 2 path for one specific Step 1 vector file instead of a whole continent batch.",
    )
    parser.add_argument(
        "--conf-factor",
        type=int,
        default=50,
        help="Cross-distance factor used in the current Step 2 run. Default: 50.",
    )
    args = parser.parse_args(argv)

    if args.vector_file:
        if args.continent is not None or args.processors is not None:
            parser.error("Use either 'continent processors' or '--vector-file', not both.")
    elif args.continent is None or args.processors is None:
        parser.error("Provide either 'continent processors' or '--vector-file'.")

    return args


def main(argv=None):
    args = parse_args(argv)
    if args.vector_file:
        process_file([args.vector_file, args.conf_factor, args.config])
    else:
        run_results_cli(
            continent_input=args.continent,
            number_of_processors=args.processors,
            config_path=args.config,
            conf_factor=args.conf_factor,
        )


if __name__ == "__main__":
    main()
