# utility functions for app.py

def check_parameter_names(config_parameters):
    '''
    ensures config parameters contain the required parameters,
    and that there are no unrecognized parameters
    min_expr, min_prop, padj_thresh, adj_method, condition,
    contrast_level, and reference_level

    returns an empty string if valid
    if invalid, returns error message
    '''

    all_parameters = {"min_expr", "min_prop", "padj_thresh", "adj_method", \
        "condition", "contrast_level", "reference_level", "use_qual_weights"}
    config_error_status = ""

    for parameter_name, parameter_value in config_parameters.items():
        if not parameter_value and parameter_name in all_parameters:
            config_error_status += \
                f"Missing value for parameter: {parameter_name}\n"
        elif parameter_name not in all_parameters:
            config_error_status += f"Unknown parameter: {parameter_name}\n"
        else:
            all_parameters.remove(parameter_name)

    # use_qual_weights is not required
    if "use_qual_weights" in all_parameters:
        all_parameters.remove("use_qual_weights")

    # if any parameters are missing, list them
    for missing_parameter in all_parameters:
        config_error_status += f"Missing parameter: {missing_parameter}\n"

    return config_error_status


def validate_parameters(config_parameters):
    '''
    returns an error message if config parameters are invalid
    otherwise, returns an empty string
    '''

    if (error_msg := check_parameter_names(config_parameters)) != "":
        return error_msg
    if type(config_parameters["min_expr"]) not in [int, float]:
        return "min_expr must be a number"
    if config_parameters["min_expr"] < 0:
        return "min_expr must be a non-negative"
    if type(config_parameters["min_prop"]) not in [int, float]:
        return "min_prop must be a number"
    if config_parameters["min_prop"] < 0:
        return "min_prop must be a non-negative"
    if type(config_parameters["padj_thresh"]) not in [int, float]:
        return "padj_thresh must be a number"
    if type(config_parameters["adj_method"]) != str:
        return "adj_method must be a string"
    if type(config_parameters["condition"]) != str:
        return "condition must be a string"
    if type(config_parameters["contrast_level"]) != str:
        return "contrast_level must be a string"
    if type(config_parameters["reference_level"]) != str:
        return "reference_level must be a string"
    if type(config_parameters["use_qual_weights"]) != bool:
        return "use_qual_weights must be a bool"

    return ""


def standardize_filename(filename):
    '''
    standardize filename to one of the following:
    counts.tsv, coldata.tsv, filter.txt or config.yml
    then return the standardized filename
    if the file name isn't recognized, return empty string

    FILE NAMING REQUIREMENTS FOR THE USER: The server recognizes filenames
        based on whether they contain one of the following unique substrings:
        "config" or ".yml": File is identified as the config.yml file
        "count": File is identified as the counts.tsv file
        "col": File is identified as the coldata.tsv file
        "filt": File is identified as the filter.txt file
    '''

    if "config" in filename or ".yml" in filename:
        filename = "config.yml"
    elif "count" in filename:
        filename = "counts.tsv"
    elif "col" in filename:
        filename = "coldata.tsv"
    elif "filt" in filename:
        filename = "filter.txt"

    return filename


def get_request_parameters(form, data_type):
    '''returns request parameters from the parameter form'''

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

