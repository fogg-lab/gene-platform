import os
import yaml
from app.job_runner import prepare_job, run_job


def update_job(directory):
    """Post an update to a job to trigger new input validation."""
    pass


def run_job(directory):
    """Start a job. Assumes input validation has already been done."""
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
