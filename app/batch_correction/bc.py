import os
import pandas as pd
import numpy as np
from IPython.display import display
import subprocess

MICROARRAY_BC_SCRIPT = "Rscript ../../rscripts/microarray_bc.r"
RNA_SEQ_BC_SCRIPT = "Rscript ../../rscripts/rnaseq_bc.r"

coldata_cols = ["sample_name", "condition", "batch"]


def call_bc(user_dir, data_type, reference_level, contrast_level):

    # Get paths to the counts and coldata files
    counts_path = f"{user_dir}counts.tsv"
    coldata_path = f"{user_dir}coldata.tsv"

    # Make sure the counts and coldata files exist
    counts_exists = os.path.isfile(counts_path)
    coldata_exists = os.path.isfile(coldata_path)
    if not counts_exists and not coldata_exists:
        return "Error: Counts and coldata files not found"
    if not counts_exists:
        return "Error: Counts file not found"
    if not coldata_exists:
        return "Error: Coldata file not found"

    # Load counts and coldata to dataframes for manipulation
    counts_cols = list(pd.read_csv(counts_path, nrows =1, sep="\t"))
    counts = pd.read_csv(counts_path, usecols = \
        [i for i in counts_cols if i.lower() != "entrez_gene_id"], sep='\t')
    coldata = pd.read_csv(coldata_path, usecols = coldata_cols, sep='\t')

    # Remove samples with levels other than the reference/contrast levels
    factor_levels = [reference_level, contrast_level]
    print(f"factor levels: {factor_levels}")
    display(coldata)
    filtered_samples = coldata.loc[coldata['condition'].isin(factor_levels)]
    display(filtered_samples)
    all_indices = set(range(0,coldata.shape[0]))
    print(f"all_indices: {all_indices}")
    samples_to_keep = set(filtered_samples.index.values)
    print(f"samples_to_keep: {samples_to_keep}")
    samples_to_drop = np.asarray(list(all_indices - samples_to_keep))
    samples_to_drop += 1    # Increment index to account for counts header
    display(samples_to_drop)
    counts.drop(counts.columns[samples_to_drop], axis=1, inplace=True)

    coldata = filtered_samples

    # Save prepared counts and coldata files for batch correction
    counts.to_csv(f"{user_dir}counts_bc-in.tsv", sep='\t', index=False)
    coldata.to_csv(f"{user_dir}coldata_bc-in.tsv", sep='\t', index=False)

    # Call the batch correction R script
    if data_type == 'microarray':
        subprocess.Popen([f"{MICROARRAY_BC_SCRIPT} {user_dir}counts_bc-in.tsv "\
                            f"{user_dir}coldata_bc-in.tsv {user_dir}"], shell=True)
    elif data_type == "RNA-Seq":
        subprocess.Popen([f"{RNA_SEQ_BC_SCRIPT} {user_dir}counts_bc-in.tsv "\
                            f"{user_dir}coldata_bc-in.tsv {user_dir}"], shell=True)

    return ""
