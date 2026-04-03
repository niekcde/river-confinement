#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse
import itertools

import numpy as np
import pandas as pd
import xarray as xr
from sklearn.decomposition import PCA

from .calc_functions import scale
from .clustering import GMM_tune, KmeansTune
from .paths import format_factor_token, load_project_paths


DEFAULT_HEIGHT_FACTORS = (2, 3, 4)
DEFAULT_CLUSTERS = tuple(range(3, 21))
DEFAULT_RANDOM_STATES = (20, 43, 50)
DEFAULT_GMM_INIT_TYPES = ("k-means++", "kmeans")
DEFAULT_GMM_MAX_ITERS = (200, 500)
DEFAULT_GMM_COV_TYPES = ("tied", "diag")
DEFAULT_GMM_INIT_COUNTS = (5, 15, 25, 50)


def _score_output_files(paths, height_factor):
    hf_token = format_factor_token(height_factor)
    return {
        "kmeans": paths.results_root / f"KmeanScores_confinement_{hf_token}.csv",
        "gmm": paths.results_root / f"GMMSCORES_confinement_{hf_token}.csv",
    }


def open_dataset_confinement_clustering(height_factor, *, cross_factor=50, config_path=None):
    paths = load_project_paths(config_path)
    cross_token = format_factor_token(cross_factor)
    hf_token = format_factor_token(height_factor)
    input_file = paths.single_smoothed_dir / f"global_{cross_token}_{hf_token}_smoothed.nc"
    if input_file.exists() is False:
        raise FileNotFoundError(
            f"Step 8 input file not found at {input_file}. "
            f"Run Step 7 smoothing for height factor {hf_token} first."
        )

    with xr.open_dataset(input_file) as ds:
        df = ds.to_dataframe().reset_index()
    return df, input_file


def prepare_confinement_clustering_dataframe(df):
    slope_cols = [
        "slope_out_normalized",
        "slope_right_normalized",
        "slope_left_normalized",
        "slope_inn_normalized",
    ]
    slope_s_cols = [
        "slope_out_normalized_smooth",
        "slope_right_normalized_smooth",
        "slope_left_normalized_smooth",
        "slope_inn_normalized_smooth",
    ]

    feature_cols = slope_cols + slope_s_cols
    df_cluster = df.dropna(subset=feature_cols).copy()
    if df_cluster.empty:
        raise ValueError(
            "Step 8 clustering input has no rows after filtering missing slope variables."
        )

    x_scaled = scale(df_cluster, feature_cols)
    n_components = min(8, x_scaled.shape[1], len(df_cluster))
    if n_components < 2:
        raise ValueError(
            "Step 8 clustering requires at least two valid PCA components after filtering."
        )

    pca = PCA(n_components=n_components)
    x_pca = pca.fit_transform(x_scaled)
    keep_components = min(3, n_components)
    pca_cols = [f"PCA{i}" for i in range(1, keep_components + 1)]
    df_cluster[pca_cols] = x_pca[:, :keep_components]
    return df_cluster, pca_cols


def _validate_clusters(clusters, n_rows):
    valid_clusters = sorted({int(c) for c in clusters if 1 < int(c) < n_rows})
    if len(valid_clusters) == 0:
        raise ValueError(
            f"No valid cluster counts remain for Step 8. "
            f"Requested clusters must be between 2 and {n_rows - 1} for the current dataset."
        )
    return valid_clusters


def _effective_sample_size(sample_size, n_rows):
    return max(2, min(int(sample_size), int(n_rows) - 1))


def run_kmeans_scores(
    df_cluster,
    cluster_cols,
    *,
    clusters,
    sample_size,
    random_states,
    cross_factor,
    height_factor,
):
    effective_sample_size = _effective_sample_size(sample_size, len(df_cluster))
    scores = KmeansTune(
        df_cluster,
        cluster_cols,
        clusters,
        effective_sample_size,
        list(random_states),
    )
    df_scores = pd.DataFrame(scores, columns=["clusters", "Silhouette", "DaviesB"])
    df_scores.insert(0, "height_factor", format_factor_token(height_factor))
    df_scores.insert(0, "cross_factor", format_factor_token(cross_factor))
    df_scores.insert(0, "n_rows", len(df_cluster))
    return df_scores


def run_gmm_scores(
    df_cluster,
    cluster_cols,
    *,
    clusters,
    sample_size,
    random_states,
    workers,
    init_types=DEFAULT_GMM_INIT_TYPES,
    max_iters=DEFAULT_GMM_MAX_ITERS,
    cov_types=DEFAULT_GMM_COV_TYPES,
    init_counts=DEFAULT_GMM_INIT_COUNTS,
    cross_factor,
    height_factor,
):
    effective_sample_size = _effective_sample_size(sample_size, len(df_cluster))
    runs = [
        (init_type, max_iter, cov_type, init_count, cluster, effective_sample_size, list(random_states))
        for init_type, max_iter, cov_type, init_count, cluster in itertools.product(
            init_types,
            max_iters,
            cov_types,
            init_counts,
            clusters,
        )
    ]

    x = df_cluster[cluster_cols]
    if workers == 1:
        scores = [GMM_tune(run, x) for run in runs]
    else:
        from joblib import Parallel, delayed

        scores = Parallel(n_jobs=workers)(delayed(GMM_tune)(run, x) for run in runs)

    df_scores = pd.DataFrame(
        scores,
        columns=[
            "init_type",
            "max_iter",
            "covType",
            "init",
            "clusters",
            "BIC",
            "AIC",
            "Silhouette",
            "DaviesB",
        ],
    )
    df_scores.insert(0, "height_factor", format_factor_token(height_factor))
    df_scores.insert(0, "cross_factor", format_factor_token(cross_factor))
    df_scores.insert(0, "n_rows", len(df_cluster))
    return df_scores


def run_confinement_clustering(
    *,
    cross_factor=50,
    height_factors=DEFAULT_HEIGHT_FACTORS,
    config_path=None,
    clusters=DEFAULT_CLUSTERS,
    sample_size=40000,
    random_states=DEFAULT_RANDOM_STATES,
    workers=1,
):
    paths = load_project_paths(config_path)
    paths.results_root.mkdir(parents=True, exist_ok=True)

    outputs = []
    for height_factor in height_factors:
        df, input_file = open_dataset_confinement_clustering(
            height_factor,
            cross_factor=cross_factor,
            config_path=config_path,
        )
        df_cluster, cluster_cols = prepare_confinement_clustering_dataframe(df)
        valid_clusters = _validate_clusters(clusters, len(df_cluster))
        output_files = _score_output_files(paths, height_factor)

        kmeans_scores = run_kmeans_scores(
            df_cluster,
            cluster_cols,
            clusters=valid_clusters,
            sample_size=sample_size,
            random_states=random_states,
            cross_factor=cross_factor,
            height_factor=height_factor,
        )
        kmeans_scores.to_csv(output_files["kmeans"], index=False)

        gmm_scores = run_gmm_scores(
            df_cluster,
            cluster_cols,
            clusters=valid_clusters,
            sample_size=sample_size,
            random_states=random_states,
            workers=workers,
            cross_factor=cross_factor,
            height_factor=height_factor,
        )
        gmm_scores.to_csv(output_files["gmm"], index=False)

        outputs.append(
            {
                "height_factor": format_factor_token(height_factor),
                "input_file": input_file,
                "kmeans_scores": output_files["kmeans"],
                "gmm_scores": output_files["gmm"],
            }
        )

    return outputs


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Step 8 entrypoint: run confinement clustering on the smoothed dataset."
    )
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--cross-factor",
        type=int,
        default=50,
        help="Cross-distance factor used in the current Step 8 run. Default: 50.",
    )
    parser.add_argument(
        "--height-factors",
        nargs="*",
        type=float,
        default=list(DEFAULT_HEIGHT_FACTORS),
        help="Height factors to cluster, for example 2 3 4. Default: 2 3 4.",
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
        help="Target sample size for silhouette and Davies-Bouldin scoring. Default: 40000.",
    )
    parser.add_argument(
        "--random-states",
        nargs="*",
        type=int,
        default=list(DEFAULT_RANDOM_STATES),
        help="Random states used during tuning. Default: 20 43 50.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of parallel workers for the GMM tuning grid. Default: 1.",
    )
    return parser.parse_args(argv)


def main_cli(argv=None):
    args = parse_args(argv)
    run_confinement_clustering(
        cross_factor=args.cross_factor,
        height_factors=args.height_factors,
        config_path=args.config,
        clusters=args.clusters,
        sample_size=args.sample_size,
        random_states=args.random_states,
        workers=args.workers,
    )


if __name__ == "__main__":
    main_cli()
