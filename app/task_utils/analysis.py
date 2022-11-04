import os.path
from flask import current_app
from redis import Redis
from rq import Queue

from app.helper import get_tsv_rows, get_analysis_confirmation_msg, StatusDict
from app.task_utils.task_runner import TaskRunner
from app.task_utils.expression_data_validation import (check_coldata_has_required_columns,
                                                       check_counts_matches_coldata,
                                                       check_factor_levels)
from app.task_utils.execute_job_async import execute_job_async


class AnalysisRunner(TaskRunner):
    """Class for preparing and running an analysis task."""

    def __init__(self, task_id, task_dir):
        super().__init__(task_id, task_dir)
        self.task_type = "analysis"
        self._input_filenames = ["counts.tsv", "coldata.tsv", "filter.txt", "config.yml"]

    MICROARRAY_SCRIPT = "dge_microarray.r"
    RNASEQ_SCRIPT = "dge_rnaseq.r"

    def execute_task(self) -> dict:
        """Queue an analysis job."""
        # get data_type param from config.yml
        data_type = self.get_config().get("data_type")

        input_dir = os.path.join(self._task_dir, "input")
        output_dir = os.path.join(self._task_dir, "output")

        log_path = os.path.join(self._task_dir, ".log")

        rscripts_path = current_app.config["RSCRIPTS_PATH"]
        if data_type == "microarray":
            script_path = os.path.join(rscripts_path, AnalysisRunner.MICROARRAY_SCRIPT)
        else:
            script_path = os.path.join(rscripts_path, AnalysisRunner.RNASEQ_SCRIPT)

        q = Queue(connection=Redis())
        cmd = f"Rscript {script_path} {input_dir} {output_dir} >>{log_path} 2>&1"
        q.enqueue(execute_job_async, args=(self._task_id, log_path, self.task_type, cmd))

        return {}

    def validate_config(self, config) -> StatusDict:
        """
        Ensures config parameters are valid
        args:
            cfg (dict): Parameters for DGE analysis.
        Returns:
            Status: Status message and errors from config validation.
        """

        status = StatusDict(status="", errors=[])

        def verify_is_numeric(key):
            value = config.get(key)
            is_numeric = False
            try:
                float(value)
                is_numeric = True
            except ValueError:
                pass
            return f"{key} must be numeric." if not is_numeric else ""

        def verify_in_range(key, min_value, max_value):
            value = config.get(key)
            in_range = min_value <= float(value) <= max_value
            return f"{key} not in range [{min_value}, {max_value}]." if not in_range else ""

        # check that all required parameters are present
        required_parameters = ("data_type", "min_expr", "min_prop", "padj_thresh", "adj_method",
                               "reference_level", "contrast_level")
        if config.get("data_type") == "microarray":
            required_parameters += ("use_qual_weights",)
        for required_param in required_parameters:
            if config.get(required_param) in (None, ""):
                status["errors"].append(f"Missing parameter: {required_param}")

        # check that the numeric parameters are valid numbers and in range
        for numeric_param in ("min_expr", "min_prop", "padj_thresh"):
            if error := verify_is_numeric(numeric_param):
                status["errors"].append(error)
            elif numeric_param in ("min_prop", "padj_thresh"):
                if error := verify_in_range(numeric_param, 0, 1):
                    status["errors"].append(error)

        # check that string parameter values are valid
        adj_methods = ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]
        if config.get("adj_method") not in adj_methods:
            status["errors"].append(f"adj_method must be one of: {adj_methods}.")
        if config.get("reference_level") == config["contrast_level"]:
            status["errors"].append("Reference level and contrast level must be different.")
        if config.get("data_type") not in ["microarray", "rnaseq"]:
            status["errors"].append("data_type must be either 'microarray' or 'rnaseq'.")

        if status["errors"]:
            return status

        status["status"] = get_analysis_confirmation_msg(config)
        return status

    def validate_task(self) -> dict:
        """Validates all input files for the task"""

        status = StatusDict(status="", errors=[])

        required_input_filenames = ("counts.tsv", "coldata.tsv", "config.yml")

        # check that all input files are present
        for filename in required_input_filenames:
            if filename not in self._input_filenames:
                status["errors"].append(f"Missing input file: {filename}")
                return status

        # check that config.yml is valid
        config = self.get_config()
        status = self.validate_config(config)

        if status["errors"]:
            return status

        # check that counts.tsv and coldata.tsv are valid
        counts = get_tsv_rows(os.path.join(self._task_dir, "counts.tsv"))
        coldata = get_tsv_rows(os.path.join(self._task_dir, "coldata.tsv"))

        reference_level, contrast_level = config["reference_level"], config["contrast_level"]

        if coldata_rows_error := check_coldata_has_required_columns(coldata):
            status["errors"].append(coldata_rows_error)
        elif coldata_levels_error := check_factor_levels(reference_level, contrast_level, coldata):
            status["errors"].append(coldata_levels_error)
        elif counts_coldata_mismatch := check_counts_matches_coldata(counts[0], coldata):
            status["errors"].append(counts_coldata_mismatch)

        return status
