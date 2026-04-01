import duckdb
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

# -------------------
# Config
# -------------------
directory = '/Volumes/PhD/confinement/results/single_values/'
src_parquet = directory + "global_50_combined.parquet"               # input parquet path string
out_parquet = directory +"areas_multiH.parquet"


H_VALUES = [5, 10, 15, 20.0, 35, 40, 45, 50]
BATCH_ROWS = 10_000                # LIMIT size (safe everywhere)
DUCKDB_THREADS = 8

# -------------------
# Core math
# -------------------
def ref_height(d, z, w):
    if d.size == 0:
        return np.nan

    mask = d <= w
    z_candidates = z[mask]

    interp = np.nan
    if d[0] <= w <= d[-1]:
        i = np.searchsorted(d, w, side="left")
        if i < d.size and d[i] == w:
            interp = z[i]
        elif 0 < i < d.size:
            d0, d1 = d[i - 1], d[i]
            z0, z1 = z[i - 1], z[i]
            if d1 != d0:
                t = (w - d0) / (d1 - d0)
                interp = z0 + t * (z1 - z0)

    if np.isfinite(interp):
        return float(min(np.min(z_candidates), interp)) if z_candidates.size else float(interp)

    return float(np.min(z_candidates)) if z_candidates.size else np.nan


# def areas_for_many_H(d, z, w, H_arr):
#     if d.size < 2:
#         return np.full(H_arr.shape, np.nan)

#     dx = np.diff(d)
#     if np.any(dx < 0):  # safety
#         idx = np.argsort(d)
#         d = d[idx]
#         z = z[idx]
#         dx = np.diff(d)

#     z_min = ref_height(d, z, w)
#     if not np.isfinite(z_min):
#         return np.full(H_arr.shape, np.nan)

#     z0, z1 = z[:-1], z[1:]
#     dx = dx.astype(np.float64)

#     out = np.empty(H_arr.shape)
#     for j, H in enumerate(H_arr):
#         z_lim = z_min + H
#         d0 = np.maximum(0.0, z_lim - z0)
#         d1 = np.maximum(0.0, z_lim - z1)
#         out[j] = np.sum((d0 + d1) * 0.5 * dx)

#     return out

import numpy as np

def areas_and_first_intersection_for_many_H(d, z, w, H_arr, *, atol=1e-12):
    if d.size < 2:
        return np.full(H_arr.shape, np.nan), np.full(H_arr.shape, np.nan)

    dx = np.diff(d)
    if np.any(dx < 0):  # safety
        idx = np.argsort(d)
        d = d[idx]
        z = z[idx]
        dx = np.diff(d)

    z_min = ref_height(d, z, w)
    if not np.isfinite(z_min):
        return np.full(H_arr.shape, np.nan), np.full(H_arr.shape, np.nan)

    z0, z1 = z[:-1], z[1:]
    x0, x1 = d[:-1], d[1:]
    dx = (x1 - x0).astype(np.float64)

    areas = np.empty(H_arr.shape, dtype=np.float64)
    x_first = np.full(H_arr.shape, np.nan, dtype=np.float64)

    for j, H in enumerate(H_arr):
        z_lim = z_min + H

        # --- area (your existing logic) ---
        d0 = np.maximum(0.0, z_lim - z0)
        d1 = np.maximum(0.0, z_lim - z1)
        areas[j] = np.sum((d0 + d1) * 0.5 * dx)

        # --- first intersection with z = z_lim ---
        s0 = z0 - z_lim
        s1 = z1 - z_lim

        # Proper crossings inside a segment
        cross = (s0 * s1 < 0.0)
        if np.any(cross):
            denom = (z1[cross] - z0[cross])
            t = (z_lim - z0[cross]) / denom
            x_cross = x0[cross] + t * (x1[cross] - x0[cross])
            x_first[j] = np.min(x_cross)
            continue

        # If no proper crossing, allow "touching" at vertices (exact endpoint on the line)
        on0 = np.isclose(s0, 0.0, atol=atol)
        on1 = np.isclose(s1, 0.0, atol=atol)
        if np.any(on0) or np.any(on1):
            candidates = []
            if np.any(on0): candidates.append(x0[on0])
            if np.any(on1): candidates.append(x1[on1])
            x_touch = np.concatenate(candidates) if candidates else np.array([], dtype=float)
            if x_touch.size:
                x_first[j] = np.min(x_touch)

    return areas, x_first



# -------------------
# Streaming via LIMIT/OFFSET
# -------------------
con = duckdb.connect()
con.execute(f"SET threads TO {DUCKDB_THREADS}")

total_rows = con.execute(
    f"SELECT COUNT(*) FROM parquet_scan('{src_parquet}')"
).fetchone()[0]

H_arr = np.asarray(H_VALUES, dtype=np.float64)
writer = None
rows_written = 0

for offset in range(0, total_rows, BATCH_ROWS):
    print(f"Processing rows {offset:,} → {min(offset+BATCH_ROWS, total_rows):,}")

    # fetch one chunk
    table = con.execute(f"""
        SELECT global_id, bendWidths, distInn, elevInn, distOut, elevOut
        FROM parquet_scan('{src_parquet}')
        LIMIT {BATCH_ROWS} OFFSET {offset}
    """).arrow()

    gids = table["global_id"].to_pylist()
    widths = table["bendWidths"].to_pylist()

    dist_inn = table["distInn"].to_pylist()
    elev_inn = table["elevInn"].to_pylist()
    dist_out = table["distOut"].to_pylist()
    elev_out = table["elevOut"].to_pylist()

    n = len(gids)
    A_inn = np.full((n, len(H_arr)), np.nan)
    A_out = np.full((n, len(H_arr)), np.nan)
    W_inn = np.full((n, len(H_arr)), np.nan)
    W_out = np.full((n, len(H_arr)), np.nan)

    for i in range(n):
        w = widths[i]
        if w is None or not np.isfinite(w) or w == 0:
            continue

        dI = np.asarray(dist_inn[i], dtype=np.float64)
        zI = np.asarray(elev_inn[i], dtype=np.float64)
        dO = np.asarray(dist_out[i], dtype=np.float64)
        zO = np.asarray(elev_out[i], dtype=np.float64)

        A_inn[i, :], W_inn[i, :] = areas_and_first_intersection_for_many_H(dI, zI, w, H_arr)
        A_out[i, :], W_out[i, :] = areas_and_first_intersection_for_many_H(dO, zO, w, H_arr)

    out_cols = {"global_id": pa.array(gids)}
    for j, H in enumerate(H_arr):
        tag = str(H).replace(".", "_")
        out_cols[f"area_inn_H{tag}"] = pa.array(A_inn[:, j])
        out_cols[f"area_out_H{tag}"] = pa.array(A_out[:, j])
        out_cols[f"W_inn_H{tag}"] = pa.array(W_inn[:, j])
        out_cols[f"W_out_H{tag}"] = pa.array(W_out[:, j])

    out_table = pa.table(out_cols)

    if writer is None:
        writer = pq.ParquetWriter(out_parquet, out_table.schema, compression="zstd")

    writer.write_table(out_table)
    rows_written += out_table.num_rows

if writer:
    writer.close()

print(f"\nDone ✔  Wrote {rows_written:,} rows → {out_parquet}")

