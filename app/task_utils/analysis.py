import os
import subprocess
from flask import current_app

from app.task_utils.task_runner import TaskRunner


class AnalysisRunner(TaskRunner):
    """Class for preparing and running an analysis task."""

    def __init__(self, task_id, task_dir):
        super().__init__(task_id, task_dir)
        self.task_type = "analysis"
        self._input_filenames = ["counts.tsv", "coldata.tsv", "filter.txt", "config.yml"]

    MICROARRAY_SCRIPT = "dge_microarray.r"
    RNASEQ_SCRIPT = "dge_rnaseq.r"

    def update_task(self):
        """Task has a new input file - perform input validation."""
        return dict(status="", warnings=[], errors=[])

    def execute_task(self):
        """Run an analysis task."""
        return dict(status="", warnings=[], errors=[])

    def validate_config(self, config):
        """
        Ensures config parameters are valid
        args:
            cfg (dict): Parameters for DGE analysis.
        Returns:
            dict: status - list of warnings under "warnings" key, errors under "errors" key
        """

        status = dict(warnings=[], errors=[])

        def in_cfg(key):
            key_in_cfg = key in config
            if not key_in_cfg:
                status["errors"].append(f"Missing parameter: {key}")
            return key_in_cfg

        def verify_is_numeric(key):
            value = config[key]
            try:
                float(value)
            except ValueError:
                status["errors"].append(f"Parameter {key} must be numeric.")

        def verify_in_range(key, min_value, max_value):
            value = config[key]
            in_range = min_value <= value <= max_value
            if not in_range:
                status["errors"].append(f"Parameter {key} must be between "
                                        f"{min_value} and {max_value}.")

        if in_cfg("min_expr"):
            verify_is_numeric("min_expr")
        if in_cfg("min_prop"):
            verify_is_numeric("min_prop")
            verify_in_range("min_prop", 0, 1)
        if in_cfg("padj_thresh"):
            verify_is_numeric("padj_thresh")
            verify_in_range("padj_thresh", 0, 1)

        adj_methods = ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]
        if in_cfg("adj_method") and config["adj_method"] not in adj_methods:
            status["errors"].append(f"Parameter adj_method must be one of {adj_methods}.")
        if config["reference_level"] == config["contrast_level"]:
            status["errors"].append("Reference level and contrast level must be different.")

        return status

    def validate_task(self) -> dict:
        """Validates all input files for the task"""
        return dict(status="", warnings=[], errors=[])

    @staticmethod
    def _call_analysis(data_type, task_dir):
        """
        Calls DGE analysis script
        Args:
            data_type (str): "microarray" or "rnaseq"
            task_dir (str): Path to task directory
        """

        input_dir = os.path.join(task_dir, "input")
        output_dir = os.path.join(task_dir, "output")

        log_path = os.path.join(task_dir, ".log")

        rscripts_path = current_app.config["RSCRIPTS_PATH"]
        if data_type == "microarray":
            script_path = os.path.join(rscripts_path, AnalysisRunner.MICROARRAY_SCRIPT)
        else:
            script_path = os.path.join(rscripts_path, AnalysisRunner.RNASEQ_SCRIPT)

        subprocess.Popen([f"{script_path} {input_dir} {output_dir} 1> {log_path} 2>& 1"],
                        shell=True)

    @staticmethod
    def _get_analysis_confirmation_msg(config_params):
        """Get analysis formula from the config parameters"""

        contrast_level = config_params["contrast_level"]
        reference_level = config_params["reference_level"]

        analysis_formula =  "<p><b>You are performing the following analysis:</b></p>\n"
        analysis_formula += f"<p><i>condition ~ (intercept) + {contrast_level}</i></p>\n\n"
        analysis_formula += f"<p>where the reference group is <i>{reference_level}</i>.\n</p>"

        return analysis_formula
