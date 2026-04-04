from __future__ import annotations

from pathlib import Path

import pandas as pd


def _parquet_engine() -> str:
    for engine in ("pyarrow", "fastparquet"):
        try:
            __import__(engine)
            return engine
        except Exception:
            continue

    raise ImportError(
        "Parquet support is not available. Install 'pyarrow' or 'fastparquet' in the active "
        "environment before running the canonical bend-table pipeline path."
    )


def read_bend_table(path: str | Path) -> pd.DataFrame:
    input_path = Path(path)
    if input_path.suffix.lower() == ".parquet":
        return pd.read_parquet(input_path, engine=_parquet_engine())
    if input_path.suffix.lower() == ".csv":
        return pd.read_csv(input_path)
    raise ValueError(f"Unsupported bend-table input format: {input_path}")


def write_bend_table(df: pd.DataFrame, path: str | Path) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_path, index=False, engine=_parquet_engine())
    return output_path
