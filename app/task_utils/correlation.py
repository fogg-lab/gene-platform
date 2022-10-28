import os
import subprocess
import pandas as pd
import numpy as np
from flask import current_app

from app.task_utils.task_runner import TaskRunner


class CorrelationRunner(TaskRunner):
    """Class for preparing and running an RNASeq correlation task."""
    def __init__(self, task_id, task_dir):
        super().__init__(task_id, task_dir)
        self.task_type = "correlation"
        self._input_filenames = ["counts.tsv", "config.yml"]

    CORR_SCRIPT = "correlation.r"

    def update_task(self):
        """Task has a new input file - perform input validation."""
        return dict(status="", warnings=[], errors=[])

    def execute_task(self):
        """Run an RNASeq correlation task"""
        return dict(status="", warnings=[], errors=[])

    def validate_config(self, config) -> dict:
        """Ensures config parameters are valid for the task"""
        return dict(status="", warnings=[], errors=[])

    def validate_task(self) -> dict:
        """Validates all input files for the task"""
        return dict(status="", warnings=[], errors=[])

    def _call_corr(self, user_dir, corr_method):
        """
        Prepare data for correlation and call R script to generate plots.
        Args:
            user_dir (string): Absolute path to directory containing input files.
            corr_method (string): Correlation method to use ('pearson', 'spearman', or 'both').
        """

        counts_path = os.path.join(user_dir, "counts.tsv")

        if not os.path.isfile(counts_path):
            return "Error: Counts file not found"

        counts_cols = list(pd.read_csv(counts_path, nrows = 1, sep="\t"))
        counts = pd.read_csv(counts_path,
            usecols = [i for i in counts_cols if i.lower() != "entrez_gene_id"], sep='\t')

        counts = counts.select_dtypes(include=np.number)

        counts.to_csv(os.path.join(user_dir, "counts_corr-in.tsv"), sep='\t', index=False)

        counts_in_path = os.path.join(user_dir, 'counts_corr-in.tsv')

        rscripts_path = current_app.config["RSCRIPTS_PATH"]
        script = os.path.join(rscripts_path, CorrelationRunner.CORR_SCRIPT)

        if corr_method not in ["pearson", "spearman", "both"]:
            raise ValueError(f"Invalid correlation method: '{corr_method}'")

        corr_methods = ["pearson", "spearman"] if corr_method == "both" else [corr_method]

        # Call R script to generate plots
        for method in corr_methods:
            cmd = [script, counts_in_path, user_dir, method]
            subprocess.Popen(cmd, cwd=user_dir, shell=True)
