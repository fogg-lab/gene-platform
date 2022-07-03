import os
import pandas as pd
import subprocess

MICROARRAY_BC_SCRIPT = "Rscript ../../rscripts/microarray_bc.r"
RNA_SEQ_BC_SCRIPT = "Rscript ../../rscripts/rnaseq_bc.r"

coldata_cols = ["sample_name", "condition", "batch"]


def call_bc(user_dir, data_type, reference_level, contrast_level):

    counts_path = f"{user_dir}counts.tsv"
    coldata_path = f"{user_dir}coldata.tsv"

    counts_exists = os.path.isfile(counts_path)
    coldata_exists = os.path.isfile(coldata_path)

    if not counts_exists and not coldata_exists:
        return "Error: Counts and coldata files not found"
    if not counts_exists:
        return "Error: Counts file not found"
    if not coldata_exists:
        return "Error: Coldata file not found"

    counts_cols = list(pd.read_csv(counts_path, nrows =1, sep="\t"))
    counts = pd.read_csv(counts_path, \
        usecols = [i for i in counts_cols if i.lower() != "entrez_gene_id"], sep='\t')
    coldata = pd.read_csv(coldata_path, usecols = coldata_cols, sep='\t')

    # TODO: remove samples with levels other than the reference/contrast levels
    # ...

    counts.to_csv(f"{user_dir}counts_bc-in.tsv", sep='\t', index=False)
    coldata.to_csv(f"{user_dir}coldata_bc-in.tsv", sep='\t', index=False)

    if data_type == 'microarray':
        subprocess.Popen([f"{MICROARRAY_BC_SCRIPT} {user_dir}counts_bc-in.tsv "\
                            f"{user_dir}coldata_bc-in.tsv {user_dir}"], shell=True)
    elif data_type == "RNA-Seq":
        subprocess.Popen([f"{RNA_SEQ_BC_SCRIPT} {user_dir}counts_bc-in.tsv "\
                            f"{user_dir}coldata_bc-in.tsv {user_dir}"], shell=True)

    return ""
