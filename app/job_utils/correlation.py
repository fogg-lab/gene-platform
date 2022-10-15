"""Class for preparing and running an RNASeq correlation job."""
import os
import subprocess
import pandas as pd
import numpy as np

from app.job_utils.job_runner import JobRunner


class CorrelationRunner(JobRunner):
    def __init__(self, job_id, job_dir):
        super().__init__(job_id, job_dir)
        self.job_type = "correlation"
        self._input_filenames = ["counts.tsv", "config.yml"]
        self._corr_script = "Rscript ../../rscripts/correlation.r"

    def update_job(self):
        """Job has a new input file - perform input validation."""
        pass


    def start_job(self):
        """Run an RNASeq correlation job"""
        pass


    def call_corr(self, user_dir, corr_method):
        """Prepare data for correlation and call R script to generate plots"""

        counts_path = os.path.join(user_dir, "counts.tsv")

        if not os.path.isfile(counts_path):
            return "Error: Counts file not found"

        counts_cols = list(pd.read_csv(counts_path, nrows = 1, sep="\t"))
        counts = pd.read_csv(counts_path,
            usecols = [i for i in counts_cols if i.lower() != "entrez_gene_id"], sep='\t')

        counts = counts.select_dtypes(include=np.number)

        counts.to_csv(os.path.join(user_dir, "counts_corr-in.tsv"), sep='\t', index=False)

        counts_in_path = os.path.join(user_dir, 'counts_corr-in.tsv')

        if corr_method != "pearson":
            subprocess.Popen([f"{self._corr_script} {counts_in_path} {user_dir} spearman"],
                             shell=True)
        if corr_method != "spearman":
            subprocess.Popen([f"{self._corr_script} {counts_in_path} {user_dir} pearson"],
                             shell=True)

        return ""
