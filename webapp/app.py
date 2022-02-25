import os
import subprocess
import shutil
import csv
import yaml
from datetime import timedelta
import tempfile
from flask import Flask, render_template, request, redirect, url_for, \
    session, Response
from flask_session.__init__ import Session


app = Flask(__name__)

app.config["SESSION_PERMANENT"] = True
app.config["SESSION_TYPE"] = "filesystem"
app.config["PERMANENT_SESSION_LIFETIME"] = timedelta(minutes=30)

Session(app)

'''
backlog:
TODO: fix params template to prefill use_qual_weights checkmark correctly
TODO: Analysis runtime counter
TODO: Add button to remove uploaded file (cancel button next to progress bar)
TODO: Parameter constraints validation
TODO: Validation:
        - Supplied column names must be checked against the values present in
            the data. If there is a mismatch (e.g., user entered “disease
            status” but the column name is “disease_condition”), the user
            should be informed of the issue and given the opportunity to
            correct the data columns and resubmit.
        - Column names for values returned by all analysis types must match
            (e.g., “log2FoldChange” (DESeq2) and “logFC” (limma) should both
            be abbreviated the same)
        - Supplied factor levels must be present in the data and user should
            be informed of errors and allowed to correct
        - Constraint: 0 < p-thresh < 1
        - Constraint: 0 < min_prop < 1
        - Constraint: do not allow navigating to parameters without uploading
            all required files
TODO: Descriptions for parameters
TODO: User story 4 (inform user of analysis to be conducted before execution)
TODO: Sort and filter
TODO: In the filter.txt error message, include the line of the error
'''

# Ensure that the current working directory is the webapp directory
# Get the path to the webapp dir from the path of this script
SCRIPT_PATH = os.path.realpath(__file__)
SCRIPT_DIR = "/".join(SCRIPT_PATH.split("/")[:-1])
os.chdir(SCRIPT_DIR)

RNA_SEQ_SCRIPT = "Rscript ../analysis/rnaseq.R"
MICROARRAY_SCRIPT = "Rscript ../analysis/microarray.R"
USER_FILES_LOCATION = "user_files"


@app.route("/")
def index():
    '''main page, first input page (file uploads)'''

    ensure_session_dir()

    user_files = os.listdir(session["user_session_dir"])
    uploads = [filename for filename in user_files if "output" not in filename]

    return render_template("uploads_form.html", uploads=uploads)


@app.route("/upload", methods=["POST"])
def upload():
    '''
    handles uploading files coldata, filter, filter and config
    one file per request
    file contents are in request.data (a bytes object)
    '''

    result = {}

    filename = request.headers.get('X_FILENAME')
    standardized_filename = standardize_filename(filename)

    if standardized_filename:
        save_temp_file(request.data, standardized_filename)
    else:
        result["error_status"] = f"Invalid filename: {filename}. \
            Filename must end with either counts.tsv, coldata.tsv, filter.txt,\
            config.yml or config.txt\n"

    if standardized_filename == "config.yml":
        result["error_status"] = check_config()

    return result


@app.route("/parameters", methods=["GET"])
def parameters():
    '''loads the parameter form'''

    parameters = parse_config()

    return render_template("parameters_form.html", params=parameters)


@app.route("/submit", methods=["POST"])
def submit():
    '''
    submit_parameters()
    submits input for analysis if valid
    generates a config file from the parameter form
    '''

    # get whether analysis is microarray or RNA-Seq
    data_type = request.form.get("data_type")

    parameters = get_request_parameters(request.form, data_type)

    generate_config(parameters)

    # remove old output
    delete_user_file("filter_output.tsv")
    delete_user_file("output.tsv")

    call_analysis(data_type)

    wait_for_output()

    return redirect(url_for("display_output"))


@app.route("/display")
def display_output():
    '''
    display output in a table
    '''

    # check here if output.tsv exists and errors.txt doesn't
    path_to_output1 =  f"{session['user_session_dir']}output.tsv"
    output1 = open(path_to_output1, encoding="UTF-8")
    reader1 = csv.reader(output1, delimiter="\t")
    output = [[elem for elem in output1] for output1 in reader1]
    output1.close()

    return render_template("results.html",output_cols=output[:1][0],output_body=output[1:])

@app.route("/filter_display")
def display_filtered():
    '''
    display filtered output in a table
    '''

    #perform the same function as displaying results, but for filtered tsv
    path_to_output2 =  f"{session['user_session_dir']}filter_output.tsv"
    output2 = open(path_to_output2, encoding="UTF-8")
    reader2 = csv.reader(output2, delimiter="\t")
    filter_output = [[elem for elem in output2] for output2 in reader2]
    output2.close()

    return render_template("filter_results.html",filter_cols=filter_output[:1][0],filter_body=filter_output[1:])


@app.route("/reset", methods=["GET"])
def reset():
    delete_user_file("counts.tsv")
    delete_user_file("coldata.tsv")
    delete_user_file("filter.txt")
    delete_user_file("config.yml")
    delete_user_file("filter_output.tsv")
    delete_user_file("output.tsv")
    return redirect(url_for("index"))


@app.context_processor
def utility_processor():
    '''defines functions that templates can use'''

    def try_float(elem):
        try:
            elem = float(elem)
        except:
            pass
        return elem

    def data_sorted_on_col(data, col_names, sort_col):
        col_index = col_names.index(sort_col)
        sorted_data = sorted(data, key=lambda line: try_float(line[col_index]))
        return sorted_data

    def data_filtered(data, col_names, filter_text, filter_col=""):
        '''
            returns rows that contain the filter text
            if filter column is unspecified, search all columns
        '''
        columns_to_search = []
        if filter_col:
            columns_to_search.append(filter_col)
        else:
            columns_to_search.append(col_names)
        filtered_data = []
        for row in data:
            row_match = False
            for colname in columns_to_search:
                col_index = col_names[colname]
                if filter_text in row[col_index]:
                    row_match = True
            if row_match:
                filtered_data.append(row)

    return dict(try_float=try_float, data_sorted_on_col=data_sorted_on_col,
        data_filtered=data_filtered)


@app.route("/getunfilteredtsv")
def get_unfiltered_tsv():
    '''download unfiltered output'''

    unfiltered_output = open(f"{session['user_session_dir']}output.tsv", \
        encoding="UTF-8")

    return Response(
        unfiltered_output,
        mimetype='text/csv',
        headers={'Content-disposition':
                'attachment; filename=output.tsv'})


@app.route("/getfilteredtsv")
def get_filtered_tsv():
    '''download filtered output'''

    filtered_output = open(f"{session['user_session_dir']}filter_output.tsv", \
        encoding="UTF-8")

    return Response(
        filtered_output,
        mimetype="text/csv",
        headers={"Content-disposition":
                "attachment; filename=filter_output.tsv"})


def save_temp_file(file_contents, filename):
    '''
    Saves uploaded file to user session directory

    Parameters:
        file_contents: bytes object
        filename: string
    '''

    # create session directory if none exists
    ensure_session_dir()

    user_file_path = f"{session['user_session_dir']}{filename}"

    if os.path.exists(user_file_path):
        os.remove(user_file_path)

    user_file = open(user_file_path, "w", encoding="UTF-8")
    lines = file_contents.split(b'\n')

    # copy the file line by line, adding a trailing newline if none exists
    for line in lines:
        user_file.write(f"{line.decode('UTF-8')}\n")

    user_file.close()


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
            config_error_status += f"Missing value for parameter: {parameter_name}\n"
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


def check_factor_levels():
    '''ensures factor levels are present in the input files'''

    pass


def call_analysis(data_type):
    '''
    calls microarray or rna-seq analysis depending on data_type
    sends the session_id as an argument
    '''

    if data_type == 'microarray':
        subprocess.Popen([f"{MICROARRAY_SCRIPT} {session['session_id']}"], \
            shell=True)
    elif data_type == 'RNA-Seq':
        subprocess.Popen([f"{RNA_SEQ_SCRIPT} {session['session_id']}"], \
                shell=True)


def read_user_file(filename):
    '''
    opens a user file for reading
    just supply the filename like "counts.tsv" for example
    returns read-only user file like counts, coldata, config, or filter
    '''
    user_file = None

    if "user_session_dir" in session:
        filepath = f"{session['user_session_dir']}{filename}"
        if os.path.isfile(filepath):
            user_file = open(filepath, "r", encoding="UTF-8")

    return user_file


def delete_user_file(filename):
    '''
    deletes a user input file if it exists
    just supply the filename like "counts.tsv" for example
    '''

    if "user_session_dir" in session:
        filepath = f"{session['user_session_dir']}{filename}"
        if os.path.isfile(filepath):
            os.remove(filepath)


def wait_for_output():
    '''busy-waits until output shows up in users session directory'''

    analysis_done = False
    while not analysis_done:
        path_to_output =  f"{session['user_session_dir']}output.tsv"

        # "2>/dev/null" suppresses the expected error output "file not found"
        output = subprocess.Popen([f"ls {path_to_output} 2>/dev/null"], \
            stdout=subprocess.PIPE, shell=True).communicate()[0]

        # results of the ls are returned in bytes, ends with newline character
        analysis_done = output == str.encode(path_to_output) + b"\n"


def cleanup_old_sessions():
    '''Clean up other sessions older than 4 hours (14400 seconds)'''

    for old_dir in os.listdir(USER_FILES_LOCATION):
        old_dir = f"{USER_FILES_LOCATION}/{old_dir}"

        if old_dir != f"{USER_FILES_LOCATION}/.gitkeep":
            get_age = "$(($(date +%s) - $(date +%s -r " + old_dir + ")))"

            # Run bash command and return the stdout output
            age_seconds = int(subprocess.Popen([f"echo {get_age}"], \
                stdout=subprocess.PIPE, shell=True).communicate()[0])

            if age_seconds > 14400:
                shutil.rmtree(old_dir)


def get_session_dir():
    if "user_session_dir" in session:
        return session["user_session_dir"]
    else:
        return False


def generate_config(config_parameters):
    '''generates a config file using user-entered parameters'''

    config_file_path = f"{session['user_session_dir']}/config.yml"
    config_file = open(config_file_path, "w", encoding="UTF-8")

    for param in config_parameters.keys():
        if param in ["adj_method", "condition", "contrast_level",\
             "reference_level"]:
            config_file.write(f"{param}: \'{config_parameters[param]}\'\n")
        else:
            config_file.write(f"{param}: {config_parameters[param]}\n")

    config_file.close()


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


def standardize_filename(filename):
    '''
    standardize filename to one of the following:
    counts.tsv, coldata.tsv, filter.txt or config.yml
    then return the standardized filename
    if the file name isn't recognized, return empty string
    '''

    if "config" in filename:
        filename = "config.yml"
    elif "count" in filename:
        filename = "counts.tsv"
    elif "col" in filename:
        filename = "coldata.tsv"
    elif "filt" in filename:
        filename = "filter.txt"

    return filename


def ensure_session_dir():
    '''
    if there is no directory for the user session, create one now
    directory path is stored in the user session variable 'user_session_dir'
    '''

    if ("user_session_dir" not in session or
        not os.path.exists(session["user_session_dir"])):
        temp_dir = tempfile.mkdtemp(dir=USER_FILES_LOCATION)
        os.chmod(temp_dir, 0o777) # give everyone rwx permission for the dir
        session["user_session_dir"] = f"{temp_dir}/"
        session["session_id"] = temp_dir.split("/")[-1:]


def parse_config():
    config_file_path = session["user_session_dir"] + "/config.yml"
    config_file = open(config_file_path)
    config_parameters = yaml.load(config_file)

    err_msg = validate_parameters(config_parameters)
    if err_msg != "":
        return err_msg
    return config_parameters


def is_numeric(string):
    num_symbols = "-.0123456789"
    return all([char in num_symbols for char in string])


def is_bool(string):
    return string in {"True", "False"}


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


def check_config():
    '''get config params in a dict then validate the params'''

    err_msg = ""

    config_params = parse_config()
    if type(config_params) == str:
        err_msg = config_params
    else:
        err_msg = validate_parameters(config_params)

    if err_msg:
        delete_user_file("config.yml")

    return err_msg


cleanup_old_sessions()
