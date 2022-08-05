import os
import pandas as pd

def combine_gdc(data_dir):
    counts_paths = []
    coldata_paths = []
    counts_dfs = []
    coldata_dfs = []

    for fname in os.listdir(data_dir):
        if "_counts.tsv" in fname:
            counts_paths.append(f"{data_dir}{fname}")
        elif "_coldata.tsv" in fname:
            coldata_paths.append(f"{data_dir}{fname}")

    for counts_path, coldata_path in zip(counts_paths, coldata_paths):
        counts_df = pd.read_csv(counts_path, sep="\t")
        coldata_df = pd.read_csv(coldata_path, sep="\t")
        counts_df.set_index("Hugo_Symbol", inplace=True)
        coldata_df.set_index("sample_name", inplace=True)
        counts_dfs.append(counts_df)
        coldata_dfs.append(coldata_df)

    combined_counts = pd.concat(counts_dfs,axis=1)
    combined_coldata = pd.concat(coldata_dfs,axis=0)

    combined_counts.to_csv(f"{data_dir}counts_preprocessed.tsv", sep="\t")
    combined_coldata.to_csv(f"{data_dir}coldata_preprocessed.tsv", sep="\t")