#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

import argparse

from .paths import load_project_paths


def build_fabdem_index(config_path: str | None = None, overwrite: bool = False) -> None:
    from .dem import build_fabdem_bounds, build_fabdem_vrt

    paths = load_project_paths(config_path)
    paths.ensure_fabdem_dirs()
    paths.validate_fabdem_inputs()

    vrt_path = build_fabdem_vrt(paths.fabdem_dir, paths.fabdem_vrt, overwrite=overwrite)
    bounds_path = build_fabdem_bounds(
        paths.fabdem_dir,
        paths.fabdem_bounds,
        dem_crs='EPSG:4326',
        overwrite=overwrite,
    )

    print(f"FABDEM directory: {paths.fabdem_dir}")
    print(f"Built FABDEM VRT: {vrt_path}")
    print(f"Built FABDEM bounds: {bounds_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Intermediate pipeline step: build the FABDEM VRT and bounds cache."
    )
    parser.add_argument(
        "--config",
        help="Path to config/paths.local.json. Defaults to config/paths.local.json when present.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Rebuild the FABDEM VRT and bounds cache even when they already exist.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_fabdem_index(config_path=args.config, overwrite=args.overwrite)
