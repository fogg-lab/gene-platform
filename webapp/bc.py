import sys
import os
import pandas as pd
import subprocess

user_dir = sys.argv[1]
data_type = sys.argv[2] # "rnaseq" or "microarray"
reference_level = sys.argv[3]
contrast_level = sys.argv[4]

MICROARRAY_BC_SCRIPT = "Rscript ../batch-correction/microarray_bc.r"
RNA_SEQ_BC_SCRIPT = "Rscript ../batch-correction/rnaseq_bc.r"

coldata_cols = ["sample_name", "condition", "batch"]

# gather all filenames for data to be batch corrected
all_files = os.listdir(user_dir)
counts_files = []
coldata_files = []

for file in all_files:
    if file == "counts.tsv":
        counts_files.append(file)
    elif file == "coldata.tsv":
        coldata_files.append(file)

if len(counts_files) != len(coldata_files):
    print("Error: number of counts files does not match number of coldata files")
    sys.exit(1)

if len(counts_files) != 1:
    print("Error: Too many input files")
    sys.exit(1)

counts_path = f"{user_dir}{counts_files[0]}"
counts_cols = list(pd.read_csv(counts_path, nrows =1, sep="\t"))
counts = pd.read_csv(counts_path, \
    usecols = [i for i in counts_cols if i.lower() != "entrez_gene_id"], sep='\t')
coldata = pd.read_csv(f"{user_dir}{coldata_files[0]}", \
    usecols = coldata_cols, sep='\t')

# TODO: remove samples with levels other than the reference/contrast levels
# ...

counts.to_csv(f"{user_dir}counts_bc-in.tsv", sep='\t', index=False)
coldata.to_csv(f"{user_dir}coldata_bc-in.tsv", sep='\t', index=False)

if data_type == 'microarray':
    subprocess.Popen([f"{MICROARRAY_BC_SCRIPT} {user_dir}counts_bc-in.tsv "\
                        f"{user_dir}coldata_bc-in.tsv {user_dir}"], shell=True)
elif data_type == "rnaseq":
    subprocess.Popen([f"{RNA_SEQ_BC_SCRIPT} {user_dir}counts_bc-in.tsv "\
                        f"{user_dir}coldata_bc-in.tsv {user_dir}"], shell=True)
