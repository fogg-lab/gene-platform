"""utility functions for app.py.  """

import os

def check_parameter_names(config_parameters):
    """
    ensures config parameters contain the required parameters,
    and that there are no unrecognized parameters
    min_expr, min_prop, padj_thresh, adj_method, condition,
    contrast_level, and reference_level

    returns an empty string if valid
    if invalid, returns error message
    """

    all_parameters = {"min_expr", "min_prop", "padj_thresh", "adj_method",
        "condition", "contrast_level", "reference_level", "use_qual_weights"}

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


def validate_parameters(config_parameters):
    """
    returns an error message if config parameters are invalid
    returns an empty string if config parameters are valid
    """

    error_msg = ""
    param_names_invalid = check_parameter_names(config_parameters)

    if param_names_invalid:
        error_msg = param_names_invalid
    else:
        if type(config_parameters["min_expr"]) not in [int, float]:
            error_msg += '"min_expr" must be a number\n'
        elif config_parameters["min_expr"] < 0:
            error_msg += '"min_expr" must be a non-negative\n'

        if type(config_parameters["min_prop"]) not in [int, float]:
            error_msg += '"min_prop" must be a number\n'
        elif not (0 <= config_parameters["min_prop"] <= 1):
            error_msg += '"min_prop" must be between 0 and 1\n'

        if type(config_parameters["padj_thresh"]) not in [int, float]:
            error_msg += '"padj_thresh" must be a number\n'
        elif not (0 <= config_parameters["padj_thresh"] <= 1):
            error_msg += '"padj_thresh" must be between 0 and 1\n'

        adj_methods = ["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]
        if not isinstance(config_parameters["adj_method"], str):
            error_msg += '"adj_method" must be a string"\n'
        elif config_parameters["adj_method"] not in adj_methods:
            error_msg += f"Unknown adjustment method: '{config_parameters['adj_method']}'\n"
            error_msg += f"Valid adj_methods: {adj_methods}\n"

        if not isinstance(config_parameters["condition"], str):
            error_msg += '"condition" must be a string"\n'

        if not isinstance(config_parameters["contrast_level"], str):
            error_msg += "contrast_level must be a string\n"

        if not isinstance(config_parameters["reference_level"], str):
            error_msg += "reference_level must be a string"

        if config_parameters["reference_level"] == \
           config_parameters["contrast_level"]:
            error_msg += \
                'Reference_level and contrast_level cannot be the same.\n'

        if "use_qual_weights" in config_parameters:
            if not isinstance(config_parameters["use_qual_weights"], bool):
                error_msg += "use_qual_weights must be either True or False.\n"

    return error_msg


def get_request_parameters(form, data_type):
    """returns request parameters from the parameter form"""

    request_parameters = {}

    # if analysis type is microarray, consider use_qual_weights
    if data_type != "RNA-Seq":
        # form.get("use_qual_weights") will initially be either 'None' or 'on'
        # it needs to be a boolean True or False
        if form.get("use_qual_weights") is None:
            request_parameters["use_qual_weights"] = False
        else:
            request_parameters["use_qual_weights"] = True

    request_parameters["min_prop"] = form.get("min_prop")
    request_parameters["min_expr"] = form.get("min_expr")
    request_parameters["adj_method"] = form.get("adj_method")
    request_parameters["condition"] = form.get("condition")
    request_parameters["contrast_level"] = form.get("contrast_level")
    request_parameters["reference_level"] = form.get("reference_level")
    request_parameters["padj_thresh"] = form.get("padj_thresh")

    return request_parameters


def check_factor_levels(config_params, coldata):
    """
    ensures factor levels are present in the input files
    if the factor levels look good, returns empty string
    """

    condition = config_params["condition"]
    contrast_level = config_params["contrast_level"]
    reference_level = config_params["reference_level"]

    condition_col_index = 0
    if condition in coldata[0]:
        condition_col_index = coldata[0].index(condition)
    else:
        return f"Condition '{condition}' not present in coldata (line 1)\n"

    contrast_level_found = False
    reference_level_found = False

    # remove header row from coldata
    coldata.pop(0)
    err_msg = ""
    for coldata_row in coldata:
        if coldata_row:
            factor_level = coldata_row[condition_col_index]
            if factor_level == contrast_level:
                contrast_level_found = True
            elif factor_level == reference_level:
                reference_level_found = True
            else:
                return f"Unknown factor level '{factor_level}'"
    if contrast_level_found is False:
        err_msg += f"Unknown contrast level '{contrast_level}'\n"
    if reference_level_found is False:
        err_msg += f"Unknown reference level '{reference_level}'\n"

    return err_msg


def check_coldata_matches_counts(counts_colnames, coldata_rows):
    """
    ensure rows in coldata match with the column names for samples in counts
    Assumes that sample names are listed on first row (header) of counts file
    Also assumes that sample names are listed in first column of coldata file,
        starting on second row of coldata file (first row after the header)
    Returns empty string if coldata and counts sample names match
    """

    samples = []

    coldata_rows.pop(0)
    for coldata_row in coldata_rows:
        if coldata_row:
            samples.append(coldata_row[0])

    if not samples:
        return "no samples are listed in the coldata file"

    # remove leading elements from counts which are not sample names
    while counts_colnames and counts_colnames[0] != samples[0]:
        counts_colnames.pop(0)

    if not counts_colnames:
        return f"sample '{samples[0]}' from coldata not found in counts file"

    # if counts and coldata match, err_msg will be empty
    err_msg = ""

    # make sure each sample name matches between counts and coldata
    while samples and counts_colnames and samples[0] == counts_colnames[0]:
        samples.pop(0)
        counts_colnames.pop(0)

    if counts_colnames and not samples:
        err_msg = f"sample '{counts_colnames[0]}' not found in coldata file"

    elif samples and not counts_colnames:
        err_msg = f"sample '{samples[0]}' not found in counts file"

    elif samples and counts_colnames and samples[0] != counts_colnames[0]:
        err_msg = f"sample '{samples[0]}' in coldata file does not match the \
            corresponding sample name '{counts_colnames[0]}' in counts file"

    return err_msg


def get_analysis_confirmation_msg(config_params):
    """get analysis formula from the config parameters"""

    condition = config_params["condition"]
    contrast_level = config_params["contrast_level"]
    reference_level = config_params["reference_level"]

    analysis_formula =  "<p><b>You are performing the following analysis:"\
                        "</b></p>\n"
    analysis_formula += f"<p><i>{condition} ~ (intercept) + "\
                        f"{contrast_level}</i></p>\n\n"
    analysis_formula += f"<p>where the reference group is <i>"\
                        f"{reference_level}</i>.\n</p>"

    return analysis_formula


def ensure_batches(coldata_rows):
    """
    ensure that all samples have a batch number
    if all samples have a batch number, returns empty string
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
                return f"Sample {sample_name} has no batch number in coldata."

    return ""


def delete_user_file(filename, user_dir):
    """
    deletes a user input file if it exists
    pass in the filename i.e "counts.tsv"
    """

    if user_dir:
        filepath = os.path.join(user_dir, filename)
        if os.path.isfile(filepath):
            os.remove(filepath)
