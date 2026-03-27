import pandas as pd
from pathlib import Path
from glob import glob
import json
import pyarrow as pa
import pyarrow.parquet as pq

# List of CSV files (can be absolute or relative paths)
directory = '/Volumes/PhD/confinement/results/single_values/'
csv_files = glob(directory + '*.csv')

def str_array_to_list(x):
    if pd.isna(x):
        return None
    if x == 99999 or x == "99999":
        return None
    if isinstance(x, list):
        return x
    try:
        return json.loads(x)
    except Exception:
        return None

# print(f"Saved Parquet file to: {output_path}")
output_file = "combined.parquet"
writer = None
count = 0
counters = {}
array_cols = ["distInn", "distOut", "elevInn", "elevOut"]

for csv in csv_files:
    
    path = Path(csv)

    # ---- derive identifiers from filename / context ----
    cont, cont_id, orthog_dist = csv.split('/')[-1].split('_')
    orthog_dist = orthog_dist.split('.')[0]  # same as split('.')[0]
    key = (cont, cont_id)
    counters.setdefault(key, 0)
    print(cont, cont_id, orthog_dist, count)
    count +=1
    # ---- stream the CSV in chunks ----
    for chunk in pd.read_csv(csv, 
                             chunksize=200_000,
                             encoding="latin1",
                            low_memory=False):
        print('\tChunk')
        # add identifier columns
        chunk["continent"] = cont
        chunk["continent_id"] = cont_id
        chunk["orthog_dist"] = orthog_dist

                # IMPORTANT: use the real column names as they appear in your CSV
        sort_cols = ["combined_reach_id", "bendDistOut"]

        # sort within the chunk (stable sort)
        chunk.sort_values(sort_cols, kind="mergesort", inplace=True, ignore_index=True)

        # unique int per row for this (continent, continent_id), continuing across chunks
        start = counters[key]
        n     = len(chunk)
        seq = pd.Series(range(start, start + n), index=chunk.index, dtype="int64")
        counters[key] = start + n

        # build a global_id (string). Safe and readable.
        chunk["global_id"] = (
            chunk["continent"].astype(str)
            + "_"
            + chunk["continent_id"].astype(str)
            + "_"
            + seq.astype(str)
        )
        
        for col in array_cols:
            chunk[col] = chunk[col].map(str_array_to_list)

        table = pa.Table.from_pandas(chunk, preserve_index=False)

        if writer is None:
            writer = pq.ParquetWriter(
                directory + f"global_{orthog_dist}_combined.parquet",
                table.schema,
                compression="snappy"
            )

        writer.write_table(table)

writer.close()