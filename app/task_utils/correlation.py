import os
import pandas as pd
import numpy as np
from redis import Redis
from rq import Queue
from flask import current_app

from app.helper import StatusDict
from app.task_utils.task_runner import TaskRunner
from app.task_utils.execute_job_async import execute_job_async


class CorrelationRunner(TaskRunner):
    """Class for preparing and running an RNASeq correlation task."""
    def __init__(self, task_id, task_dir):
        super().__init__(task_id, task_dir)
        self.task_type = "correlation"
        self._input_filenames = ["counts.tsv", "config.yml"]

    CORR_SCRIPT = "correlation.r"

    def execute_task(self):
        """Run an RNASeq correlation task"""
        cfg = self.get_config()
        corr_method = cfg.get("corr_method")
        return self._call_corr(corr_method)

    def validate_config(self, config) -> dict:
        """Ensures config parameters are valid for the task"""
        return StatusDict(status="", errors=[])

    def validate_task(self) -> dict:
        """Validates all input files for the task"""
        return StatusDict(status="", errors=[])

    def _call_corr(self, corr_method):
        """
        Prepare data for correlation and call R script to generate plots.
        Args:
            user_dir (string): Absolute path to directory containing input files.
            corr_method (string): Correlation method to use ('pearson', 'spearman', or 'both').
        """

        in_dir = os.path.join(self._task_dir, "input")
        out_dir = os.path.join(self._task_dir, "output")
        counts_path = os.path.join(in_dir, "counts.tsv")

        if not os.path.isfile(counts_path):
            err_msg = "Error: Counts file not found"
            return StatusDict(status=err_msg, errors=[err_msg])

        counts_cols = list(pd.read_csv(counts_path, nrows = 1, sep="\t"))
        counts = pd.read_csv(counts_path,
            usecols = [i for i in counts_cols if i.lower() != "entrez_gene_id"], sep='\t')

        counts = counts.select_dtypes(include=np.number)

        counts_in_path = os.path.join(in_dir, 'counts_corr-in.tsv')
        counts.to_csv(counts_in_path, sep='\t', index=False)

        scripts_path = current_app.config["SCRIPTS_PATH"]
        script = os.path.join(scripts_path, CorrelationRunner.CORR_SCRIPT)

        if corr_method not in ["pearson", "spearman", "both"]:
            raise ValueError(f"Invalid correlation method: '{corr_method}'")

        corr_methods = ["pearson", "spearman"] if corr_method == "both" else [corr_method]

        cmd = [f"Rscript {script} {counts_in_path} {out_dir} {method}" for method in corr_methods]
        cmd = " && ".join(cmd)
        cmd = f"bash -c '{cmd}'"

        log_path = os.path.join(self._task_dir, ".log")

        q = Queue(connection=Redis())
        q.enqueue(execute_job_async, args=(self._task_id, log_path, "correlation", cmd))

