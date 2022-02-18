import sys

import os
import subprocess
import shutil
import csv
from datetime import timedelta
import tempfile
from werkzeug.utils import secure_filename
from flask import Flask, render_template, request, redirect, url_for, \
    session, Response
from flask_session.__init__ import Session
from flask_dropzone import Dropzone


app = Flask(__name__)
dropzone = Dropzone(app)

app.config["SESSION_PERMANENT"] = True
app.config["SESSION_TYPE"] = "filesystem"
app.config["PERMANENT_SESSION_LIFETIME"] = timedelta(minutes=30)

Session(app)


#TODO: nav buttons, file validation


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

    files = os.listdir(USER_FILES_LOCATION)
    return render_template("uploads_form.html", files=files)


@app.route("/uploadcounts", methods=["POST"])
def upload_counts():
    '''handles uploading a counts file'''

    uploaded_file = request.files["file"]
    filename = secure_filename(uploaded_file.filename)
    file_ext = os.path.splitext(filename)[1]

    if file_ext != ".tsv":
        return "Invalid file extension", 400

    save_temp_file(uploaded_file, "counts.tsv")

    return redirect(url_for("index"))


@app.route("/uploadcoldata", methods=["POST"])
def upload_coldata():
    '''handles uploading a coldata file'''

    uploaded_file = request.files["file"]
    filename = secure_filename(uploaded_file.filename)
    file_ext = os.path.splitext(filename)[1]

    if file_ext != ".tsv":
        return "Invalid file extension", 400

    save_temp_file(uploaded_file, "coldata.tsv")

    return redirect(url_for("index"))


@app.route("/uploadfilter", methods=["POST"])
def upload_filter():
    '''handles uploading a filter file'''

    uploaded_file = request.files["file"]
    filename = secure_filename(uploaded_file.filename)
    file_ext = os.path.splitext(filename)[1]

    if file_ext != ".txt":
        return "Invalid file extension", 400

    save_temp_file(uploaded_file, "filter.txt")

    return redirect(url_for("index"))


@app.route("/uploadconfig", methods=["POST"])
def upload_config():
    '''handles uploading a config file'''

    uploaded_file = request.files["file"]
    filename = secure_filename(uploaded_file.filename)
    file_ext = os.path.splitext(filename)[1]

    if file_ext not in {".txt", ".yml"}:
        return "Invalid file extension", 400

    save_temp_file(uploaded_file, "config.yml")

    config_validity = validate_config()
    if config_validity != "valid":
        return config_validity, 400

    return redirect(url_for("index"))


@app.route("/parameters", methods=["GET"])
def parameters_form():
    '''loads the parameter form'''

    parameters = get_config_parameters()

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

    if not data_type:
        pass    # error - blank datatype. todo: inform user

    parameters = get_request_parameters(request.form, data_type)

    generate_config(parameters)

    call_analysis(data_type)
    remove_old_output()

    wait_for_output()

    return redirect(url_for("display_output"))


@app.route("/display")
def display_output():
    '''
    display output in a table
    '''

    # check here if output.tsv exists and errors.txt doesn"t
    path_to_output =  f"{session['user_session_dir']}output.tsv"
    output = open(path_to_output, encoding="UTF-8")
    reader = csv.reader(output, delimiter="\t")
    rows = [[elem for elem in row] for row in reader]
    output.close()

    return render_template("results.html", col_titles=rows[:1], info=rows[1:])


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


def save_temp_file(file, filename):
    '''
    Takes a user"s file and copies it into a temp directory on the server
    directory path is stored in the user session variable 'user_session_dir'
    '''

    if ("user_session_dir" not in session or
            not os.path.exists(session["user_session_dir"])):
        temp_dir = tempfile.mkdtemp(dir=USER_FILES_LOCATION)
        os.chmod(temp_dir, 0o777) # give everyone rwx permission for the dir
        session["user_session_dir"] = f"{temp_dir}/"
        session["session_id"] = temp_dir.split("/")[-1:]

    user_file_path = f"{session['user_session_dir']}{filename}"

    if os.path.exists(user_file_path):
        os.remove(user_file_path)

    user_file = open(user_file_path, "w", encoding="UTF-8")

    # copy the file line by line, adding a trailing newline if none exists
    lines = file.readlines()

    for line in lines:
        user_file.write(line.decode("UTF-8"))

    user_file.close()


def validate_config():
    '''
    validates config.yml, returns the status, deletes the file if invalid
    
    ensures config has all needed parameters:
    min_expr, min_prop, padj_thresh, adj_method, condition, contrast_level,
    and reference_level

    returns either "valid", "Missing value for parameter: xxxx", or
    "Missing parameter: xxxx" where "xxxx" is replaced w/ the name of the param
    '''

    config_parameters = get_config_parameters()

    all_parameters = {"min_expr", "min_prop", "padj_thresh", "adj_method", \
        "condition", "contrast_level", "reference_level", "use_qual_weights"}
    config_status = ""

    for parameter_name, parameter_value in config_parameters.items():
        if not parameter_value and parameter_name in all_parameters:
            config_status += f"Missing value for parameter: {parameter_name}\n"
        elif parameter_name not in all_parameters:
            config_status += f"Unknown parameter: {parameter_name}\n"
        else:
            all_parameters.remove(parameter_name)

    # use_qual_weights is not required
    if "use_qual_weights" in all_parameters:
        all_parameters.remove("use_qual_weights")

    # if any parameters are missing, list them
    for missing_parameter in all_parameters:
        config_status += f"Missing parameter: {missing_parameter}\n"

    if not config_status:
        config_status = "valid"
    else:
        delete_user_file("config.yml")

    return config_status


def validate_filter():
    '''ensures that filter has one gene per line'''

    filter_file = read_user_file("filter.txt")
    filter_status = "valid"

    for line in filter_file:
        word_list = line.split()
        if len(word_list) > 1:
            filter_status = "invalid: not one gene per line"

    if filter_file:
        filter_file.close()

    return filter_status


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
    user_file = None

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
    '''Clean up other sessions older than one day (86400 seconds)'''

    for old_dir in os.listdir(USER_FILES_LOCATION):
        old_dir = f"{USER_FILES_LOCATION}/{old_dir}"

        if old_dir != f"{USER_FILES_LOCATION}/.gitkeep":
            get_age = "$(($(date +%s) - $(date +%s -r " + old_dir + ")))"

            # Run bash command and return the stdout output
            age_seconds = int(subprocess.Popen([f"echo {get_age}"], \
                stdout=subprocess.PIPE, shell=True).communicate()[0])

            if age_seconds > 28800:
                shutil.rmtree(old_dir)


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


def get_config_parameters():
    '''parses config to return a dictionary of parameters with set values'''
    config_file = None

    if "user_session_dir" in session:
        print('Hello world!', file=sys.stderr)
        config_filepath = f"{session['user_session_dir']}config.yml"
        print(config_filepath, file=sys.stderr)
        if os.path.isfile(config_filepath):
            print('Hello world2!', file=sys.stderr)
            config_file = open(config_filepath, "r", encoding="UTF-8")

    config_parameters = {}

    if config_file:
        config_lines = config_file.readlines()
        for line in config_lines:
            # remove quotes and newline characters
            line = line.translate(str.maketrans("", "", "\n\'\""))
            if line:
                parameter_name, parameter_value = tuple(line.split(": "))
                config_parameters[parameter_name] = parameter_value

    return config_parameters


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


def remove_old_output():
    '''removes old output files from the current session'''
    filtered_output_path = f"{session['user_session_dir']}filter_output.tsv"
    unfiltered_output_path = f"{session['user_session_dir']}output.tsv"
    filtered_output_exists = os.path.exists(filtered_output_path)
    unfiltered_output_exists = os.path.exists(unfiltered_output_path)
    if filtered_output_exists:
        os.remove(filtered_output_path)
    if unfiltered_output_exists:
        os.remove(unfiltered_output_path)


cleanup_old_sessions()
