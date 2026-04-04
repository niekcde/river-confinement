#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse

from .clustering_confinement import (
    DEFAULT_CLUSTERS,
    DEFAULT_HEIGHT_FACTORS,
    DEFAULT_RANDOM_STATES,
    run_confinement_clustering,
)
from .spatial_smoothing import run_spatial_smoothing


def run_step7_8_confinement_workflow(
    *,
    cross_factor=50,
    height_factors=DEFAULT_HEIGHT_FACTORS,
    config_path=None,
    smoothing_workers=6,
    clustering_workers=1,
    continents=None,
    clusters=DEFAULT_CLUSTERS,
    sample_size=40000,
    random_states=DEFAULT_RANDOM_STATES,
):
    smoothing_outputs = []
    for height_factor in height_factors:
        smoothing_outputs.append(
            run_spatial_smoothing(
                cross_factor=cross_factor,
                height_factor=height_factor,
                config_path=config_path,
                workers=smoothing_workers,
                continents=continents,
            )
        )

    clustering_outputs = run_confinement_clustering(
        cross_factor=cross_factor,
        height_factors=height_factors,
        config_path=config_path,
        clusters=clusters,
        sample_size=sample_size,
        random_states=random_states,
        workers=clustering_workers,
    )

    return {
        "smoothing_outputs": smoothing_outputs,
        "clustering_outputs": clustering_outputs,
    }


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Step 7-8 workflow: smooth and cluster the required confinement height factors together."
    )
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--cross-factor",
        type=int,
        default=50,
        help="Cross-distance factor used in the current Step 7-8 workflow. Default: 50.",
    )
    parser.add_argument(
        "--height-factors",
        nargs="*",
        type=float,
        default=list(DEFAULT_HEIGHT_FACTORS),
        help="Height factors to smooth and cluster, for example 2 3 4. Default: 2 3 4.",
    )
    parser.add_argument(
        "--continents",
        nargs="*",
        help="Optional subset of continent codes to smooth before clustering, for example oc af.",
    )
    parser.add_argument(
        "--smoothing-workers",
        type=int,
        default=6,
        help="Number of continent-level worker processes for Step 7 smoothing. Default: 6.",
    )
    parser.add_argument(
        "--clustering-workers",
        type=int,
        default=1,
        help="Number of parallel workers for Step 8 clustering. Default: 1.",
    )
    parser.add_argument(
        "--clusters",
        nargs="*",
        type=int,
        default=list(DEFAULT_CLUSTERS),
        help="Cluster counts to evaluate. Default: 3 through 20.",
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=40000,
        help="Target sample size for Step 8 scoring. Default: 40000.",
    )
    parser.add_argument(
        "--random-states",
        nargs="*",
        type=int,
        default=list(DEFAULT_RANDOM_STATES),
        help="Random states used during Step 8 tuning. Default: 20 43 50.",
    )
    return parser.parse_args(argv)


def main_cli(argv=None):
    args = parse_args(argv)
    run_step7_8_confinement_workflow(
        cross_factor=args.cross_factor,
        height_factors=args.height_factors,
        config_path=args.config,
        smoothing_workers=args.smoothing_workers,
        clustering_workers=args.clustering_workers,
        continents=args.continents,
        clusters=args.clusters,
        sample_size=args.sample_size,
        random_states=args.random_states,
    )


if __name__ == "__main__":
    main_cli()
