import os
import subprocess
from pathlib import Path
import pandas as pd
import numpy as np
from redis import Redis
from rq import Queue
from flask import current_app

from app.helper import get_tsv_rows, StatusDict
from app.task_utils.task_runner import TaskRunner
from app.task_utils.execute_job_async import execute_job_async


class BatchCorrectionRunner(TaskRunner):
    """Class for preparing and running a batch correction task."""
    def __init__(self, task_id, task_dir):
        super().__init__(task_id, task_dir)
        self.task_type = "batch"
        self._input_filenames = ["counts.tsv", "coldata.tsv", "config.yml"]

    BC_SCRIPT = "batch_correction.r"

    def execute_task(self) -> StatusDict:
        """Queue a batch correction job."""
        input_dir = Path(self._task_dir) / "input"
        cfg = self.get_config()
        data_type = cfg.get("data_type")
        reference_level = cfg.get("reference_level")
        contrast_level = cfg.get("contrast_level")
        return BatchCorrectionRunner._call_batch_correction(self._task_id, input_dir, data_type,
                                                            reference_level, contrast_level)

    def validate_config(self, config) -> StatusDict:
        """Ensures config parameters are valid for the task"""
        return StatusDict(status="", errors=[])

    def validate_task(self) -> StatusDict:
        """Validates all input files for the task"""
        return StatusDict(status="", errors=[])

    @staticmethod
    def _call_batch_correction(task_id, directory, data_type, reference_level,
                               contrast_level) -> StatusDict:
        """
        Calls the R script to run batch correction.
        Args:
            task_id (str): The task id.
            directory (str): Absolute path to directory containing input files.
            data_type (str): Type of expression data, either "microarray" or "rnaseq".
            reference_level (str): Reference level for the samples, i.e "normal".
            contrast_level (str): Contrast level for the samples, i.e "tumor".
        Returns:
            StatusDict: A status message and errors, if any.
        """

        status = StatusDict(status="", errors=[])

        # Get paths to the counts and coldata files
        counts_path = os.path.join(directory, "counts.tsv")
        coldata_path = os.path.join(directory, "coldata.tsv")

        # Read in counts and coldata dataframes
        counts_cols = list(pd.read_csv(counts_path, nrows=1, sep="\t"))
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

        counts_in = os.path.join(directory, "counts_bc-in.tsv")
        coldata_in = os.path.join(directory, "coldata_bc-in.tsv")

        # Save prepared counts and coldata files for batch correction
        counts.to_csv(counts_in, sep='\t', index=False)
        coldata.to_csv(coldata_in, sep='\t', index=False)

        # Call the batch correction R script
        script = os.path.join(current_app.config["SCRIPTS_PATH"], BatchCorrectionRunner.BC_SCRIPT)
        log_path = os.path.join(directory, ".log")
        cmd = f"{script} {counts_in} {coldata_in} {directory} {data_type}"

        if not status.get("errors"):
            q = Queue(connection=Redis())
            q.enqueue(execute_job_async, args=(task_id, log_path, "batch correction", cmd))

        return status

    def check_bc_input_files(self):
        """
        Ensures batch correction input files are present.
        Args:
            directory (string): Absolute path to directory containing input files.
        Returns:
            string: Empty string if input files are present, error message otherwise.
        """

        input_dir = Path(self._task_dir) / "input"

        counts_path = os.path.join(input_dir, "counts.tsv")
        coldata_path = os.path.join(input_dir, "coldata.tsv")

        counts_exists = os.path.isfile(counts_path)
        coldata_exists = os.path.isfile(coldata_path)

        if not counts_exists and not coldata_exists:
            return "Error: Counts and coldata files not found"
        if not counts_exists:
            return "Error: Counts file not found"
        if not coldata_exists:
            return "Error: Coldata file not found"

        return ""

    def check_batch_correction_coldata(self):
        """
        Ensures coldata has batches.
        Returns:
            string: Empty string if coldata has batches, error message otherwise.
        """

        input_dir = os.path.join(self._task_dir, "input")
        coldata_path = os.path.join(input_dir, "coldata.tsv")
        coldata_rows = get_tsv_rows(coldata_path)
        status_msg = BatchCorrectionRunner._ensure_batches(coldata_rows)

        if status_msg:
            os.remove(coldata_path)

        return status_msg

    @staticmethod
    def _ensure_batches(coldata_rows):
        """
        Ensures all samples have a batch number.
        Args:
            coldata_rows (list): List of sample rows from a coldata file.
        Returns:
            string: Empty string if all samples have a batch number, error message otherwise.
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
                    return f"Sample '{sample_name}' has no batch number in coldata."

        return ""
