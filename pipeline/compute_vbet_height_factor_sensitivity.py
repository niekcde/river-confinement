"""Recalculate confinement ratios for VBET-matched manuscript bends.

This is deliberately a sample-scale calculation: it reads the original
cross-sectional profiles only for bends in the VBET match table, obtains their
stored ``conFactor`` values from the factor-2 NetCDF, and recalculates the
one-sided ratios for selected height factors.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import xarray as xr

from pipeline.calc_functions import confinement_values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matches",
        default="results/vbet_validation/pnw_dgo_bend_matches.parquet",
        type=Path,
    )
    parser.add_argument(
        "--source",
        default="/Volumes/PhD/confinement/results/single_values/global_50_combined.parquet",
        type=Path,
    )
    parser.add_argument(
        "--factor-2-netcdf",
        default="/Volumes/PhD/confinement/results/single_values/global_50_02_conf.nc",
        type=Path,
    )
    parser.add_argument(
        "--output",
        default="results/vbet_validation/pnw_bend_height_factor_sensitivity.parquet",
        type=Path,
    )
    parser.add_argument("--match-distance-m", type=float, default=150)
    return parser.parse_args()


def ratio_for_factor(row: pd.Series, height_factor: float) -> tuple[float, float]:
    """Return outer and inner ratios, matching the production calculation."""
    try:
        _, _, _, _, _, _, er_out, er_inn = confinement_values(
            row.elevOut,
            row.elevInn,
            row.distOut,
            row.distInn,
            row.bendWidths,
            row.bendMaxWidths,
            row.conFactor * height_factor,
        )
    except (TypeError, ValueError, IndexError):
        return np.nan, np.nan
    return er_out, er_inn


def main() -> None:
    args = parse_args()
    for path in (args.matches, args.source, args.factor_2_netcdf):
        if not path.exists():
            raise FileNotFoundError(path)

    matches = pd.read_parquet(args.matches)
    matches = matches.loc[matches.match_distance_m <= args.match_distance_m].copy()
    bends = (
        matches.groupby("bendID", as_index=False)
        .agg(
            combined_reach_id=("combined_reach_id_x", "first"),
            bendDistOut=("bendDistOut", "first"),
            n_vbet_dgos=("dgoid", "size"),
            vbet_integrated_width_m=("integrated_width", "median"),
            vbet_channel_width_m=("vbet_channel_width", "median"),
            paper_ER_out_stored=("ER_out", "first"),
            paper_ER_inn_stored=("ER_inn", "first"),
        )
    )

    # Avoid materialising North America's raw-profile columns: join only the
    # requested bend keys before projecting the profile arrays.
    with duckdb.connect() as con:
        con.register("selected_bends", bends[["combined_reach_id", "bendDistOut"]])
        source = con.execute(
            f'''SELECT s."Unnamed: 0" AS source_row, s.global_id,
                       s.combined_reach_id, s.bendDistOut, s.bendWidths,
                       s.bendMaxWidths, s.distOut, s.distInn, s.elevOut, s.elevInn,
                       s.n_chan_max, s.n_chan_mod, s.widthRatio, s.wse, s.slope,
                       s.facc, s.strm_order, s.networkGroup
                FROM read_parquet('{args.source}') AS s
                INNER JOIN selected_bends AS b
                    ON s.combined_reach_id = b.combined_reach_id
                   AND s.bendDistOut = b.bendDistOut
            '''
        ).fetchdf()

    # conFactor was retained in the factor-2 confinement output. The original
    # row-number field is only unique within input tiles, so use the stable
    # bend key instead.
    with xr.open_dataset(args.factor_2_netcdf) as ds:
        factors = ds[["combined_reach_id", "bendDistOut", "conFactor"]].to_dataframe().reset_index(drop=True)
    factors = factors.drop_duplicates(["combined_reach_id", "bendDistOut"])
    source = source.merge(
        factors,
        on=["combined_reach_id", "bendDistOut"],
        how="left",
        validate="one_to_one",
    )
    result = bends.merge(
        source,
        on=["combined_reach_id", "bendDistOut"],
        how="left",
        validate="one_to_one",
    )

    for factor in (1.5, 2.0, 3.0):
        values = result.apply(ratio_for_factor, axis=1, height_factor=factor)
        result[f"ER_out_h{factor:g}"] = [value[0] for value in values]
        result[f"ER_inn_h{factor:g}"] = [value[1] for value in values]
        result[f"paper_ratio_h{factor:g}"] = (
            result[f"ER_out_h{factor:g}"] + result[f"ER_inn_h{factor:g}"]
        )
        result[f"paper_valley_width_h{factor:g}_m"] = (
            result[f"paper_ratio_h{factor:g}"] * result.bendWidths
        )

    result["vbet_ratio"] = result.vbet_integrated_width_m / result.bendWidths
    result["paper_ratio_stored"] = result.paper_ER_out_stored + result.paper_ER_inn_stored
    result["factor2_ratio_difference"] = result.paper_ratio_h2 - result.paper_ratio_stored

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(args.output, index=False)
    valid = result.factor2_ratio_difference.dropna()
    print(f"Wrote {len(result):,} bends to {args.output}")
    print(f"Source profiles recovered: {result.source_row.notna().sum():,}/{len(result):,}")
    print(
        "Factor-2 reproduction (recomputed minus stored ratio): "
        f"median={valid.median():.3g}; median absolute difference={valid.abs().median():.3g}; n={len(valid):,}"
    )


if __name__ == "__main__":
    main()
