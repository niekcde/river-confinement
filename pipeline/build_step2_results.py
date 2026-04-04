#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse

from .step2 import run_results_cli


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 2 entrypoint: build orthogonals, sample DEM profiles, and write results/all outputs."
    )
    parser.add_argument("continent", help="Continent code to process, for example 'af' or 'eu'.")
    parser.add_argument("processors", type=int, help="Number of worker processes to use.")
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--conf-factor",
        type=int,
        default=50,
        help="Cross-distance factor used in the current Step 2 run. Default: 50.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    run_results_cli(
        continent_input=args.continent,
        number_of_processors=args.processors,
        config_path=args.config,
        conf_factor=args.conf_factor,
    )


if __name__ == "__main__":
    main()
