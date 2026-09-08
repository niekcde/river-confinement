"""Directional Step 7 smoothing; no network-wide shortest-path tables."""
from dataclasses import dataclass

import numpy as np
import pandas as pd


ATTRIBUTES = [f"slope_{side}_normalized" for side in ("left", "right", "out", "inn")]
OUTPUT_ATTRIBUTES = [f"{a}_smooth" for a in ATTRIBUTES] + [f"{a}_smoothSTD" for a in ATTRIBUTES]


def prepare_bends(df):
    df = df.copy()
    group = ["file", "networkGraph", "combined_reach_id"]
    if df[group + ["bendDistOut"]].isna().any().any():
        raise ValueError("Missing bend ordering or network identifiers")
    if df.duplicated(group + ["bendDistOut"]).any():
        raise ValueError("Tied bendDistOut values: bend ordering is ambiguous")
    df["bendRank"] = df.groupby(group)["bendDistOut"].rank(ascending=False).astype(int)
    df["bendID"] = df["combined_reach_id"].astype(int).astype(str) + "_" + df["bendRank"].astype(str)
    return df.sort_values(group + ["bendRank"])


@dataclass
class BendTopology:
    lengths: np.ndarray
    upstream: np.ndarray
    downstream: np.ndarray
    paths: np.ndarray
    diagnostics: dict


def build_topology(df, max_neighbors=4):
    """Rows must have ranks from prepare_bends; paths contain row positions.

    Upstream follows the single combined_reach_up reference already selected
    upstream of Step 7. We do not traverse all tributaries at a confluence.
    """
    if max_neighbors < 1:
        raise ValueError("max_neighbors must be positive")
    lengths = df["bendLen"].to_numpy(dtype=float)
    if not np.all(np.isfinite(lengths) & (lengths > 0)):
        raise ValueError("bendLen must be finite and strictly positive")
    work = df.reset_index(drop=True)
    keys = ["file", "networkGraph", "combined_reach_id"]
    reaches = {key: rows.sort_values("bendRank").index.to_numpy()
               for key, rows in work.groupby(keys, sort=False)}
    up = np.full(len(df), -1, dtype=int)
    dn = up.copy()
    missing = 0
    for key, positions in reaches.items():
        for col in ("combined_reach_up", "combined_reach_dn"):
            if work.loc[positions, col].nunique(dropna=False) > 1:
                raise ValueError(f"Inconsistent {col} within reach {key}")
        up[positions[1:]] = positions[:-1]
        dn[positions[:-1]] = positions[1:]
        for pos, col, target, endpoint in (
            (positions[0], "combined_reach_up", up, -1),
            (positions[-1], "combined_reach_dn", dn, 0),
        ):
            ref = work.at[pos, col]
            if pd.notna(ref):
                neighbor = reaches.get((key[0], key[1], ref))
                if neighbor is None:
                    missing += 1
                else:
                    target[pos] = neighbor[endpoint]
    paths = np.full((len(df), 2, max_neighbors), -1, dtype=int)
    repeats = 0
    for focal in range(len(df)):
        seen = {focal}
        for direction, adjacency in enumerate((up, dn)):
            current = focal
            for hop in range(max_neighbors):
                candidate = adjacency[current]
                if candidate < 0:
                    break
                if candidate in seen:
                    repeats += 1
                    break
                seen.add(candidate)
                paths[focal, direction, hop] = candidate
                current = candidate
    if repeats:
        raise ValueError(f"Detected {repeats} repeated-node walks; inspect network connectivity")
    return BendTopology(lengths, up, dn, paths, {"missing_or_outside_network_references": missing})


def candidate_weights(topology, focal, neighbors=3, alpha=0.75, length_floor=True):
    if not isinstance(neighbors, (int, np.integer)) or not 1 <= neighbors <= topology.paths.shape[2]:
        raise ValueError("neighbors must be a positive integer within cached path depth")
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be finite and positive")
    length = topology.lengths[focal]
    indices, actual, effective, directions = [focal], [0.0], [0.0], ["self"]
    for direction, name in enumerate(("up", "down")):
        previous = length
        previous_effective = length
        distance = distance_effective = 0.0
        for candidate in topology.paths[focal, direction, :neighbors]:
            if candidate < 0:
                break
            current = topology.lengths[candidate]
            current_effective = max(current, length) if length_floor else current
            distance += (previous + current) / 2
            distance_effective += (previous_effective + current_effective) / 2
            indices.append(candidate)
            actual.append(distance)
            effective.append(distance_effective)
            directions.append(name)
            previous, previous_effective = current, current_effective
    weights = np.exp(-0.5 * (np.asarray(effective) / (alpha * length)) ** 2)
    weights /= weights.sum()
    return np.asarray(indices), np.asarray(actual), np.asarray(effective), weights, directions


def smooth_local(df, topology=None, *, neighbors=3, alpha=0.75, length_floor=True):
    """Preserve legacy NaN propagation and sample-count STD correction."""
    if topology is None:
        topology = build_topology(df, neighbors)
    values = df[ATTRIBUTES].to_numpy(dtype=float)
    if np.isinf(values).any():
        raise ValueError("Infinite normalized slope values; inspect Step 7 normalization")
    result = np.empty((len(df), 8))
    diagnostics = []
    for focal in range(len(df)):
        ids, actual, effective, weights, directions = candidate_weights(
            topology, focal, neighbors, alpha, length_floor)
        mean = np.average(values[ids], weights=weights, axis=0)
        std = (np.zeros(4) if len(ids) == 1 else
               np.sqrt((weights @ (values[ids] - mean) ** 2) / ((len(ids) - 1) / len(ids))))
        result[focal] = np.concatenate((mean, std))
        diagnostics.append((directions.count("up"), directions.count("down"),
                            actual.max(), effective.max(), 1 - weights[0]))
    output = df.copy()
    output[OUTPUT_ATTRIBUTES] = result
    output[["smoothing_n_up", "smoothing_n_down", "smoothing_max_actual_distance",
            "smoothing_max_effective_distance", "smoothing_neighbor_weight"]] = diagnostics
    return output
