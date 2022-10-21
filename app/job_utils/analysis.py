import os
import subprocess
from flask import current_app

from app.job_utils.job_runner import JobRunner


class AnalysisRunner(JobRunner):
    """Class for preparing and running an analysis job."""

    def __init__(self, job_id, job_dir):
        super().__init__(job_id, job_dir)
        self.job_type = "analysis"
        self._input_filenames = ["counts.tsv", "coldata.tsv", "filter.txt", "config.yml"]

    MICROARRAY_SCRIPT = "dge_microarray.r"
    RNASEQ_SCRIPT = "dge_rnaseq.r"

    def update_job(self):
        """Job has a new input file - perform input validation."""
        pass

    def start_job(self):
        """Run an analysis job."""
        pass

    @staticmethod
    def validate_config(cfg):
        """
        Ensures config parameters are valid
        args:
            cfg (dict): Parameters for DGE analysis.
        Returns:
            dict: status - list of warnings under "warnings" key, errors under "errors" key
        """

        status = dict(warnings=[], errors=[])

        def in_cfg(key):
            key_in_cfg = key in cfg
            if not key_in_cfg:
                status["errors"].append(f"Missing parameter: {key}")
            return key_in_cfg

        def is_numeric(key):
            value = cfg[key]
            try:
                float(value)
                return True
            except ValueError:
                status["errors"].append(f"Parameter {key} must be numeric.")
                return False

        def in_range(key, min_value, max_value):
            value = cfg[key]
            in_range = min_value <= value <= max_value
            if not in_range:
                status["errors"].append(f"Parameter {key} must be between "
                                        f"{min_value} and {max_value}.")
            return in_range

        if in_cfg("min_expr"):
            is_numeric("min_expr")
        if in_cfg("min_prop"):
            is_numeric("min_prop")
            in_range("min_prop", 0, 1)
        if in_cfg("padj_thresh"):
            is_numeric("padj_thresh")
            in_range("padj_thresh", 0, 1)

        adj_methods = ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]
        if in_cfg("adj_method") and cfg["adj_method"] not in adj_methods:
            status["errors"].append(f"Parameter adj_method must be one of {adj_methods}.")
        if cfg["reference_level"] == cfg["contrast_level"]:
            status["errors"].append("Reference level and contrast level must be different.")

        return status

    @staticmethod
    def call_analysis(data_type, job_dir):
        """
        Calls DGE analysis script
        Args:
            data_type (str): "microarray" or "rnaseq"
            job_dir (str): Path to job directory
        """

        input_dir = os.path.join(job_dir, "input")
        output_dir = os.path.join(job_dir, "output")

        log_path = os.path.join(job_dir, ".log")

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
