import os
from redis import Redis
from rq import Queue
from flask import current_app

from app.helper import StatusDict
from app.task_utils.task_runner import TaskRunner
from app.task_utils.execute_job_async import execute_job_async

class NormalizationRunner(TaskRunner):
    """Class for preparing and running a normalization task."""
    def __init__(self, task_id, task_dir):
        super().__init__(task_id, task_dir)
        self.task_type = "normalization"
        self._input_filenames = ["counts.tsv", "coldata.tsv", "config.yml"]

    NORM_SCRIPT = "normalize.r"

    def execute_task(self):
        """Run a normalization task"""
        task_dir = self._task_dir
        cfg = self.get_config()
        method = cfg.get("method")
        return NormalizationRunner._call_normalization(self._task_id, task_dir, method)

    def validate_config(self, config) -> dict:
        """Ensures config parameters are valid for the task"""
        if "method" not in config:
            return StatusDict(status="Invalid config", errors=["Missing method"])
        elif config["method"] not in ["tmm", "mrn"]:
            return StatusDict(status="Invalid config", errors=["Invalid method"])
        return StatusDict(status="", errors=[])

    def validate_task(self) -> dict:
        """Validates all input files for the task"""
        return StatusDict(status="", errors=[])

    @staticmethod
    def _call_normalization(task_id, task_dir, method) -> StatusDict:
        """
        Calls the R script to run normalization.
        Args:
            task_id (str): The task id.
            task_dir (str): Absolute path to task directory.
            method (str): The normalization method to use.
        Returns:
            StatusDict: A status message and errors, if any.
        """
        status = StatusDict(status="", errors=[])
        in_dir = os.path.join(task_dir, "input")
        out_dir = os.path.join(task_dir, "output")

        # Make sure counts and coldata files
        counts_exists = os.path.isfile(os.path.join(in_dir, "counts.tsv"))
        coldata_exists = os.path.isfile(os.path.join(in_dir, "coldata.tsv"))

        err_msg = ""
        if not counts_exists and not coldata_exists:
            err_msg = "Error: Counts and coldata files not found"
        elif not counts_exists:
            err_msg = "Error: Counts file not found"
        elif not coldata_exists:
            err_msg = "Error: Coldata file not found"
        if err_msg:
            return StatusDict(status=err_msg, errors=[err_msg])

        script = os.path.join(current_app.config["SCRIPTS_PATH"], NormalizationRunner.NORM_SCRIPT)
        log_path = os.path.join(task_dir, ".log")
        cmd = f"Rscript {script} {in_dir} {out_dir} {method}"

        q = Queue(connection=Redis())
        q.enqueue(execute_job_async, args=(task_id, log_path, "normalization", cmd))

        return status
