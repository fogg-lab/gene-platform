"""
Functions for preparing and running analysis.
Used by the job runner module.
"""
import os
import subprocess
import csv
import time
import yaml
from flask import current_app

INPUT_FNAMES = ["counts.tsv", "coldata.tsv", "filter.txt", "config.yml"]
MICROARRAY_SCRIPT = "dge_microarray.r"
RNASEQ_SCRIPT = "dge_rnaseq.r"

def update_job(directory):
    """Job has a new input file - perform input validation."""
    pass


def start_job(directory):
    """Run an analysis job."""
    pass


def get_analysis_confirmation_msg(config_params):
    """get analysis formula from the config parameters"""

    contrast_level = config_params["contrast_level"]
    reference_level = config_params["reference_level"]

    analysis_formula =  "<p><b>You are performing the following analysis:</b></p>\n"
    analysis_formula += f"<p><i>condition ~ (intercept) + {contrast_level}</i></p>\n\n"
    analysis_formula += f"<p>where the reference group is <i>{reference_level}</i>.\n</p>"

    return analysis_formula


def check_analysis_config(config_path):
    """
    Validates analysis config parameters
    Args:
        path (string): Absolute path to config file.
    Returns:
        str: Error message if invalid, empty string if valid
    """

    config_params = dict()

    if os.path.exists(config_path):
        config_file = open(config_path, encoding="UTF-8")
        config_params = yaml.load(config_file, Loader=yaml.FullLoader)
        config_file.close()

    status_msg = check_dge_analysis_parameter_names(config_params)
    if len(status_msg) == 0:
        status_msg = check_dge_analysis_parameters(config_params)

    if len(status_msg) > 0:
        os.remove(config_path)

    return status_msg


def check_dge_analysis_parameter_names(config_parameters):
    """
    Ensures config parameters contain the required parameters,
    and that there are no unrecognized parameters

    Returns an empty string if valid
    Returns error message if invalid
    """

    all_parameters = {"min_expr", "min_prop", "padj_thresh", "adj_method",
                      "contrast_level", "reference_level", "use_qual_weights"}

    config_error_status = ""
    unknown_params = []
    params_missing_value = []
    missing_params = []

    for parameter_name, parameter_value in config_parameters.items():
        if parameter_name in all_parameters and parameter_value in [None, ""]:
            params_missing_value.append(parameter_name)

        if parameter_name not in all_parameters:
            unknown_params.append(parameter_name)
        else:
            all_parameters.remove(parameter_name)

    # if any parameters are missing, list them
    if all_parameters is not None:
        for missing_parameter in all_parameters:
            if missing_parameter != "use_qual_weights":
                missing_params.append(missing_parameter)

    for param in params_missing_value:
        config_error_status += f"Missing value for parameter: {param}\n"
    for param in unknown_params:
        config_error_status += f"Unknown parameter: {param}\n"
    for param in missing_params:
        config_error_status += f"Missing parameter: {param}\n"

    return config_error_status


def check_dge_analysis_parameters(config_parameters):
    """
    Ensures config parameters are valid
    args:
        config_parameters (dict): Parameters for DGE analysis.
    Returns:
        string: Error message if invalid, empty string if valid
    """

    status_msg = ""

    if type(config_parameters["min_expr"]) not in [int, float]:
        status_msg += '"min_expr" must be a number\n'
    elif config_parameters["min_expr"] < 0:
        status_msg += '"min_expr" must be a non-negative\n'

    if type(config_parameters["min_prop"]) not in [int, float]:
        status_msg += '"min_prop" must be a number\n'
    elif not 0 <= config_parameters["min_prop"] <= 1:
        status_msg += '"min_prop" must be between 0 and 1\n'

    if type(config_parameters["padj_thresh"]) not in [int, float]:
        status_msg += '"padj_thresh" must be a number\n'
    elif not 0 <= config_parameters["padj_thresh"] <= 1:
        status_msg += '"padj_thresh" must be between 0 and 1\n'

    adj_methods = ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]
    if not isinstance(config_parameters["adj_method"], str):
        status_msg += '"adj_method" must be a string"\n'
    elif config_parameters["adj_method"] not in adj_methods:
        status_msg += f"Unknown adjustment method: '{config_parameters['adj_method']}'\n"
        status_msg += f"Valid adj_methods: {adj_methods}\n"

    if not isinstance(config_parameters["contrast_level"], str):
        status_msg += "contrast_level must be a string\n"

    if not isinstance(config_parameters["reference_level"], str):
        status_msg += "reference_level must be a string"

    if config_parameters["reference_level"] == config_parameters["contrast_level"]:
        status_msg += 'Reference_level and contrast_level cannot be the same.\n'

    if "use_qual_weights" in config_parameters:
        if not isinstance(config_parameters["use_qual_weights"], bool):
            status_msg += "use_qual_weights must be either True or False.\n"

    return status_msg


def call_analysis(data_type, job_dir):
    """
    calls microarray or rnaseq analysis depending on data type
    sends the session path as an argument, redirects output to a log file
    """

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