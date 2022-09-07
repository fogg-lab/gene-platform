"""
Functions for preparing and running batch correction.
Used by the job runner module.
"""
import os
import subprocess
import pandas as pd
import numpy as np
import csv
import time
import yaml
from flask import current_app
from app import helper

INPUT_FNAMES = ["counts.tsv", "coldata.tsv", "config.yml"]


def update_job(directory):
    """Job has a new input file - perform input validation."""
    pass


def start_job(directory):
    """"Run a batch correction job"""
    pass


def call_batch_correction(directory, data_type, reference_level, contrast_level):
    """
    Calls the R script to run batch correction.
    Args:
        directory (string): Absolute path to directory containing input files.
        data_type (string): Type of expression data, either "microarray" or "rnaseq".
        reference_level (string): Reference level for the samples, i.e "normal".
        contrast_level (string): Contrast level for the samples, i.e "tumor".
    """

    # Get paths to the counts and coldata files
    counts_path = os.path.join(directory, "counts.tsv")
    coldata_path = os.path.join(directory, "coldata.tsv")

    # Read in counts and coldata dataframes
    counts_cols = list(pd.read_csv(counts_path, nrows =1, sep="\t"))
    counts = pd.read_csv(counts_path, usecols =
        [i for i in counts_cols if i.lower() != "entrez_gene_id"], sep='\t')
    counts.rename(columns=lambda x: "symbol" if "symbol" in x.lower() else x, inplace=True)

    coldata_cols = ["sample_name", "condition", "batch"]
    coldata = pd.read_csv(coldata_path, usecols = coldata_cols, sep='\t')
    coldata.columns = coldata.columns.str.lower()

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

    counts_in_path = os.path.join(directory, "counts_bc-in.tsv")
    coldata_in_path = os.path.join(directory, "coldata_bc-in.tsv")

    # Save prepared counts and coldata files for batch correction
    counts.to_csv(counts_in_path, sep='\t', index=False)
    coldata.to_csv(coldata_in_path, sep='\t', index=False)

    # Call the batch correction R script
    bc_script = os.path.join(current_app.config["RSCRIPTS_PATH"], "batch_correction.r")
    subprocess.Popen(
        [f"{bc_script} {counts_in_path} {coldata_in_path} {directory} {data_type}"],
        shell=True
    )


def check_bc_input_files(directory):
    """
    Ensures batch correction input files are present.
    Args:
        directory (string): Absolute path to directory containing input files.
    Returns:
        string: Empty string if input files are present, error message otherwise.
    """

    counts_path = os.path.join(directory, "counts.tsv")
    coldata_path = os.path.join(directory, "coldata.tsv")

    counts_exists = os.path.isfile(counts_path)
    coldata_exists = os.path.isfile(coldata_path)

    if not counts_exists and not coldata_exists:
        return "Error: Counts and coldata files not found"
    if not counts_exists:
        return "Error: Counts file not found"
    if not coldata_exists:
        return "Error: Coldata file not found"

    return ""


def check_batch_correction_coldata(directory):
    """
    Ensures coldata has batches
    Args:
        directory: string
    Returns:
        string
    """

    coldata_path = os.path.join(directory, "coldata.tsv")
    coldata = helper.get_tsv_rows(coldata_path)
    status_msg = ensure_batches(coldata)

    if status_msg:
        os.remove(coldata_path)

    return status_msg


def ensure_batches(coldata_rows):
    """
    ensure that all samples have a batch number
    if all samples have a batch number, returns empty string
    """

    # check header row to make sure there is a batch column
    header = coldata_rows.pop(0)
    if "batch" not in header:
        return "Batch column not found in coldata file"
    batch_col_index = header.index("batch")

    # check if all samples have a batch number
    for coldata_row in coldata_rows:
        if coldata_row:
            if not coldata_row[batch_col_index].isnumeric():
                sample_name = coldata_row[0]
                return f"Sample {sample_name} has no batch number in coldata."

    return ""
