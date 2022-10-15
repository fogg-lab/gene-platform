"""Class for preparing and running an analysis job."""
import os
import subprocess
import time
from flask import current_app

from app.job_utils.job_runner import JobRunner


class AnalysisRunner(JobRunner):
    def __init__(self, job_id, job_dir):
        super().__init__(job_id, job_dir)
        self.job_type = "analysis"
        self._input_filenames = ["counts.tsv", "coldata.tsv", "filter.txt", "config.yml"]


    def update_job(self):
        """Job has a new input file - perform input validation."""
        pass


    def start_job(self):
        """Run an analysis job."""
        pass


MICROARRAY_SCRIPT = "dge_microarray.r"
RNASEQ_SCRIPT = "dge_rnaseq.r"

def get_analysis_confirmation_msg(config_params):
    """Get analysis formula from the config parameters"""

    contrast_level = config_params["contrast_level"]
    reference_level = config_params["reference_level"]

    analysis_formula =  "<p><b>You are performing the following analysis:</b></p>\n"
    analysis_formula += f"<p><i>condition ~ (intercept) + {contrast_level}</i></p>\n\n"
    analysis_formula += f"<p>where the reference group is <i>{reference_level}</i>.\n</p>"

    return analysis_formula


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


def call_analysis(data_type, job_dir):
    """calls microarray or rnaseq analysis depending on data type"""

    input_dir = os.path.join(job_dir, "input")
    output_dir = os.path.join(job_dir, "output")

    log_path = os.path.join(job_dir, ".log")

    rscripts_path = current_app.config["RSCRIPTS_PATH"]
    if data_type == "microarray":
        script_dir = os.path.join(rscripts_path, MICROARRAY_SCRIPT)
    else:
        script_dir = os.path.join(rscripts_path, RNASEQ_SCRIPT)

    subprocess.Popen([f"{script_dir} {input_dir} {output_dir} 1> {log_path} 2>& 1"],
                     shell=True)


def wait_for_output():
    """waits until output shows up in users session directory"""

    user_dir = common.Job.get_dir(job_id)

    unfiltered_output_path = os.path.join(user_dir, "output.tsv")
    filt_output_path = os.path.join(user_dir, "filter_output.tsv")
    filter_path = os.path.join(user_dir, "filter.txt")

    is_output = False
    while not is_output:
        is_output = os.path.exists(unfiltered_output_path)
        if is_output and os.path.exists(filter_path):
            is_output = os.path.exists(filt_output_path)

        if not is_output:
            time.sleep(0.5)

""" OLD CODE in submit():
# remove old output
for old_output_file in ["filter_output.tsv", "output.tsv"]:
    if os.path.isfile(old_output_file):
        os.remove(old_output_file)

for filename in os.listdir(session["user_session_dir"]):
    if filename.endswith(".png"):
        helpers.delete_user_file(filename, common.Job.get_dir(job_id))

call_analysis(data_type)

wait_for_output()
"""