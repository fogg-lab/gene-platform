import os
import pandas as pd
import numpy as np
import subprocess

BC_SCRIPT = "Rscript ../../rscripts/batch_correction.r"

coldata_cols = ["sample_name", "condition", "batch"]

def start_job(job_dir):
    pass

def check_bc_input_files(user_dir):
    """Make sure counts and coldata files exist"""

    counts_path = os.path.join(user_dir, "counts.tsv")
    coldata_path = os.path.join(user_dir, "coldata.tsv")

    counts_exists = os.path.isfile(counts_path)
    coldata_exists = os.path.isfile(coldata_path)

    if not counts_exists and not coldata_exists:
        return "Error: Counts and coldata files not found"
    if not counts_exists:
        return "Error: Counts file not found"
    if not coldata_exists:
        return "Error: Coldata file not found"

    return ""


def call_bc(user_dir, data_type, reference_level, contrast_level):

    # Get paths to the counts and coldata files
    counts_path = os.path.join(user_dir, "counts.tsv")
    coldata_path = os.path.join(user_dir, "coldata.tsv")

    # Read in counts and coldata dataframes
    counts_cols = list(pd.read_csv(counts_path, nrows =1, sep="\t"))
    counts = pd.read_csv(counts_path, usecols =
        [i for i in counts_cols if i.lower() != "entrez_gene_id"], sep='\t')
    counts.rename(columns=lambda x: "symbol" if "symbol" in x.lower() else x, inplace=True)

    coldata = pd.read_csv(coldata_path, usecols = coldata_cols, sep='\t')
    coldata.columns= coldata.columns.str.lower()

    # Remove samples with levels other than the reference/contrast levels
    factor_levels = [reference_level, contrast_level]
    filtered_samples = coldata.loc[coldata['condition'].isin(factor_levels)]
    all_indices = set(range(0,coldata.shape[0]))
    samples_to_keep = set(filtered_samples.index.values)
    samples_to_drop = np.asarray(list(all_indices - samples_to_keep))
    samples_to_drop += 1    # Increment index to account for counts header
    print(samples_to_drop)
    if len(samples_to_drop) > 0:
        counts.drop(counts.columns[samples_to_drop], axis=1, inplace=True)

    coldata = filtered_samples

    counts_in_path = os.path.join(user_dir, "counts_bc-in.tsv")
    coldata_in_path = os.path.join(user_dir, "coldata_bc-in.tsv")

    # Save prepared counts and coldata files for batch correction
    counts.to_csv(counts_in_path, sep='\t', index=False)
    coldata.to_csv(coldata_in_path, sep='\t', index=False)

    # Call the batch correction R script
    subprocess.Popen(
        [f"{BC_SCRIPT} {counts_in_path} {coldata_in_path} {user_dir} {data_type}"],
        shell=True
    )


def check_bc_coldata():
    """Ensures coldata has batches"""

    coldata = common.get_tsv_rows("coldata.tsv")
    err_msg = helpers.ensure_batches(coldata)

    if err_msg:
        helpers.delete_user_file("coldata.tsv",  common.get_session_dir())

    return err_msg
