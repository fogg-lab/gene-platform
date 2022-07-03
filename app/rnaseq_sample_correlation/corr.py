import os
import pandas as pd
import subprocess

CORR_SCRIPT = "Rscript ../../rscripts/correlation.r"

def call_corr(user_dir, corr_method):

    counts_path = f"{user_dir}counts.tsv"

    if not os.path.isfile(counts_path):
        return "Error: Counts file not found"

    counts_cols = list(pd.read_csv(counts_path, nrows = 1, sep="\t"))
    counts = pd.read_csv(counts_path, \
        usecols = [i for i in counts_cols if i.lower() != "entrez_gene_id"], sep='\t')

    counts = counts._get_numeric_data()

    counts.to_csv(f"{user_dir}counts_corr-in.tsv", sep='\t', index=False)

    if corr_method != "pearson":
        subprocess.Popen(
            [f"{CORR_SCRIPT} {user_dir}counts_corr-in.tsv {user_dir} spearman"],
             shell=True)
    if corr_method != "spearman":
        subprocess.Popen(
            [f"{CORR_SCRIPT} {user_dir}counts_corr-in.tsv {user_dir} pearson"],
             shell=True)

    return ""
