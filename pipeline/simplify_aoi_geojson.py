"""Simplify a single-feature GeoJSON AOI to a target upload size.

The script decimates the extremely dense source boundary at a regular interval.
It is intended for an upload outline, not for measurement or analysis. Keep the
original official boundary alongside the simplified copy.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def _decimate_ring(ring: list[list[float]], step: int) -> list[list[float]]:
    if ring[0] != ring[-1]:
        raise ValueError("Polygon ring must be closed.")
    simplified = ring[:-1][::step]
    if len(simplified) < 3:
        raise ValueError("Simplification would leave fewer than three ring vertices.")
    return simplified + [simplified[0]]


def _clip_ring_to_max_latitude(ring: list[list[float]], max_latitude: float) -> list[list[float]]:
    """Clip a closed longitude/latitude ring to the half-plane south of a latitude."""
    if ring[0] != ring[-1]:
        raise ValueError("Polygon ring must be closed.")
    output: list[list[float]] = []
    for start, end in zip(ring[:-1], ring[1:]):
        start_inside, end_inside = start[1] <= max_latitude, end[1] <= max_latitude
        if start_inside:
            output.append(start)
        if start_inside != end_inside:
            fraction = (max_latitude - start[1]) / (end[1] - start[1])
            output.append([start[0] + fraction * (end[0] - start[0]), max_latitude])
    if len(output) < 3:
        raise ValueError("Latitude clip removed the polygon.")
    return output + [output[0]]


def simplify_to_size(
    source: Path, destination: Path, max_bytes: int, max_latitude: float | None = None
) -> tuple[int, int]:
    with source.open(encoding="utf-8") as handle:
        data = json.load(handle)
    features = data.get("features", [])
    if len(features) != 1 or features[0].get("geometry", {}).get("type") != "Polygon":
        raise ValueError("This helper currently expects one Polygon feature.")
    original_rings = features[0]["geometry"]["coordinates"]
    if max_latitude is not None:
        original_rings = [
            _clip_ring_to_max_latitude(ring, max_latitude) for ring in original_rings
        ]
    step = max(1, int(source.stat().st_size / max_bytes) + 1)
    while True:
        data["features"][0]["geometry"]["coordinates"] = [
            _decimate_ring(ring, step) for ring in original_rings
        ]
        rendered = json.dumps(data, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
        if len(rendered) <= max_bytes:
            break
        step += 1
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(rendered)
    written = json.loads(rendered)
    vertex_count = sum(
        len(ring) for ring in written["features"][0]["geometry"]["coordinates"]
    )
    return vertex_count, step


def main() -> None:
    parser = argparse.ArgumentParser(description="Simplify a Polygon GeoJSON AOI for an upload-size limit.")
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--max-mb", type=float, default=4.75, help="Target file size below the uploader limit. Default: 4.75 MB")
    parser.add_argument(
        "--max-latitude",
        type=float,
        help="Optionally clip the AOI south of this WGS84 latitude before simplifying.",
    )
    args = parser.parse_args()
    vertices, step = simplify_to_size(
        args.source, args.destination, int(args.max_mb * 1024 * 1024), args.max_latitude
    )
    print(f"{args.destination}: {vertices} vertices; retained every {step}th source vertex")


if __name__ == "__main__":
    main()
