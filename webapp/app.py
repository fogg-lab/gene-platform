import os
import subprocess
import shutil
import csv
from datetime import timedelta
import tempfile
import yaml
from webapp import helpers
from flask import Flask, render_template, request, redirect, url_for, \
    session, Response, jsonify, send_from_directory
from flask_session.__init__ import Session

app = Flask(__name__)

app.config["SESSION_PERMANENT"] = True
app.config["SESSION_TYPE"] = "filesystem"
app.config["PERMANENT_SESSION_LIFETIME"] = timedelta(minutes=30)

Session(app)

'''
backlog:
TODO: Descriptions for the columns included in the output\
TODO: Validate string parameters
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
    '''main page'''
    return render_template("home.html", title="Welcome!")
    
@app.route("/uploadsetup")
def uploadsetup():
    '''first input page (file uploads)'''

    ensure_session_dir()

    user_files = os.listdir(session["user_session_dir"])
    uploads = {}
    for filename in ["counts.tsv", "coldata.tsv", "filter.txt", "config.yml"]:
        uploads[filename] = False
    for filename in user_files:
        if filename in uploads:
            uploads[filename] = True

    return render_template("uploads_form.html", uploads=uploads, title="Uploads")


@app.route("/upload", methods=["POST"])
def upload():
    '''
    handles uploading counts, coldata, filter, and config file
    one file per request
    file contents are in request.data (a bytes object)
    '''

    result = {}

    filename = request.headers.get('X_FILENAME')

    save_temp_file(request.data, filename)
    
    if filename == "config.yml":
        result["error_status"] = check_config()

    return jsonify(result)


@app.route("/cancelupload", methods=["POST"])
def cancelupload():
    '''
    removes uploaded file from user session directory
    '''

    filename = request.form.get("filename")
    delete_user_file(filename)

    return f"{filename} upload cancelled"


@app.route("/parameters", methods=["GET"])
def parameters():
    '''loads the parameter form'''
    desc = []
    with open('static/files/description.txt') as infile:
        line = infile.readline()
        desc.append(line)
        while line:
            line = infile.readline()
            desc.append(line)

    return render_template("parameters_form.html", params=parse_config(), title="Parameters", 
        padj_thresh = desc[0], min_prop = desc[1], min_expr = desc[2], adj_method = desc[3], condition = desc[4],
        contrast_level = desc[5], reference_level = desc[6])


@app.route("/submit", methods=["POST"])
def submit():
    '''
    submit_parameters()
    submits input for analysis if valid
    generates a config file from the parameter form
    '''

    # get whether analysis is microarray or RNA-Seq
    data_type = request.form.get("data_type")

    # remove old output
    delete_user_file("filter_output.tsv")
    delete_user_file("output.tsv")

    call_analysis(data_type)

    wait_for_output()

    return "analysis completed"


@app.route("/confirmsubmission", methods=["POST"])
def confirm_submission():
    '''validate input and display formula before submission'''

    # get whether analysis is microarray or RNA-Seq, and get params
    data_type = request.form.get("data_type")
    params = helpers.get_request_parameters(request.form, data_type)

    # generate config file from the form parameters
    generate_config(params)

    # validate the config, counts and coldata
    config_file_error = check_config()

    counts_colnames = get_tsv_rows("counts.tsv")[0]

    coldata_counts_match_error = helpers.check_coldata_matches_counts(
        counts_colnames, get_tsv_rows("coldata.tsv"))
    
    factor_levels_error = helpers.check_factor_levels(
        params, get_tsv_rows("coldata.tsv"))

    # get the analysis formula to display for the user
    confirmation_message = ""

    if coldata_counts_match_error:
        confirmation_message += f"<p>Error: {coldata_counts_match_error}</p>"
    if factor_levels_error:
        confirmation_message += f"<p>Error: {factor_levels_error}</p>"
    if config_file_error:
        confirmation_message += f"<p>Error: {config_file_error}</p>"

    if not confirmation_message:
        confirmation_message = helpers.get_confirmation_message(params)

    return confirmation_message


@app.route("/getconsoleoutput")
def get_console_output():
    '''
    returns the contents of the log file in user session directory
    the log file contains terminal output from the analysis script
    '''

    log = read_user_file("log")

    if not log:
        return ("", 204)

    return Response(log, mimetype='text/plain')


@app.route("/display")
def display_output():
    '''
    display output in a table
    '''

    # check here if output.tsv exists and errors.txt doesn't
    unfiltered_output_path =  f"{session['user_session_dir']}output.tsv"
    unfiltered_output_file = open(unfiltered_output_path, encoding="UTF-8")
    output_reader = csv.reader(unfiltered_output_file, delimiter="\t")
    unfiltered_output = list(output_reader)
    unfiltered_output_file.close()

    filter_output_path = f"{session['user_session_dir']}{'filter_output.tsv'}"
    filter_output_exists = os.path.exists(filter_output_path)
    filtered_data = None
    if filter_output_exists:
        filter_output_file = open(filter_output_path, encoding="UTF-8")
        output_reader = csv.reader(filter_output_file, delimiter="\t")
        filtered_output = list(output_reader)
        filter_output_file.close()
        filtered_data = filtered_output[1:]

    return render_template("results.html", cols = unfiltered_output[:1][0],\
         data=unfiltered_output[1:], filtered_data=filtered_data)


@app.route("/plots")
def plots():
    '''
    display plots
    '''

    session_dir = f"{session['user_session_dir']}"
    plot_filenames = []

    if session_dir:
        for filename in os.listdir(session_dir):
            if ".png" in filename:
                plot_filenames.append(filename)
    
    return render_template("plots.html", plot_filenames=plot_filenames)


@app.route('/getplot/<filename>')
def getplot(filename):
    return send_from_directory(session['user_session_dir'], filename)


@app.route("/reset", methods=["GET"])
def reset():
    '''Deletes files from a users session'''

    delete_user_file("counts.tsv")
    delete_user_file("coldata.tsv")
    delete_user_file("filter.txt")
    delete_user_file("config.yml")
    delete_user_file("filter_output.tsv")
    delete_user_file("output.tsv")

    return redirect(url_for("index"))


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

    filtered_output = open(f"{session['user_session_dir']}filter_output.tsv",\
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


def call_analysis(data_type):
    '''
    calls microarray or rna-seq analysis depending on data_type
    sends the session_id as an argument
    '''

    log = get_session_dir() + "log"

    if data_type == 'microarray':
        subprocess.Popen([f"{MICROARRAY_SCRIPT} {session['session_id']} "\
                            f"1> {log} 2>& 1"], shell=True)
    elif data_type == 'RNA-Seq':
        subprocess.Popen([f"{RNA_SEQ_SCRIPT} {session['session_id']} "\
                            f"1> {log} 2>& 1"], shell=True)


def read_user_file(filename):
    '''
    opens a user file for reading
    pass in the filename i.e "counts.tsv"
    returns lines from a file in the users session directory
    '''
    user_file = None

    session_dir = get_session_dir()
    if session_dir:
        filepath = f"{session_dir}{filename}"
        if os.path.isfile(filepath):
            user_file = open(filepath, "r", encoding="UTF-8")

    return user_file


def delete_user_file(filename):
    '''
    deletes a user input file if it exists
    pass in the filename i.e "counts.tsv"
    '''

    session_dir = get_session_dir()
    if session_dir:
        filepath = f"{session_dir}{filename}"
        if os.path.isfile(filepath):
            os.remove(filepath)


def wait_for_output():
    '''busy-waits until output shows up in users session directory'''

    analysis_done = False
    while not analysis_done:
        unfilt_output_path =  f"{session['user_session_dir']}output.tsv"
        filt_output_path =  f"{session['user_session_dir']}filter_output.tsv"
        filter_path = f"{session['user_session_dir']}filter.txt"
        is_output = os.path.exists(unfilt_output_path)
        if is_output and os.path.exists(filter_path):
            is_output = os.path.exists(filt_output_path)

        # results of the ls are returned in bytes, ends with newline character
        analysis_done = is_output


def cleanup_old_sessions():
    '''Clean up other sessions older than 4 hours (14400 seconds)'''

    for old_dir in os.listdir(USER_FILES_LOCATION):
        old_dir = f"{USER_FILES_LOCATION}/{old_dir}"

        if old_dir != f"{USER_FILES_LOCATION}/.gitkeep":
            get_age = "$(($(date +%s) - $(date +%s -r " + old_dir + ")))"

            # Run bash command and return the stdout output
            age_seconds = subprocess.Popen([f"echo {get_age}"], \
                stdout=subprocess.PIPE, shell=True).communicate()[0]

            try:
                age_seconds = int(age_seconds)
                if age_seconds > 14400:
                    shutil.rmtree(old_dir)
            except ValueError:
                pass


def get_session_dir():
    ''' returns session dir if it exists, otherwise returns False '''

    if "user_session_dir" in session:
        return session["user_session_dir"]
    else:
        return ""


def generate_config(config_parameters):
    '''generates a config file using user-entered parameters'''

    config_file_path = f"{session['user_session_dir']}config.yml"
    config_file = open(config_file_path, "w", encoding="UTF-8")

    for param in config_parameters.keys():
        if param in ["adj_method", "condition", "contrast_level",\
             "reference_level"]:
            config_file.write(f"{param}: \'{config_parameters[param]}\'\n")
        else:
            config_file.write(f"{param}: {config_parameters[param]}\n")

    config_file.close()


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
    ''' parses config into a dict of parameters '''

    session_dir = get_session_dir()
    config_file_path = session_dir + "config.yml"
    config_params = {}

    if session_dir and os.path.exists(config_file_path):
        config_file = open(config_file_path, encoding="UTF-8")
        config_params = yaml.safe_load(config_file)
        config_file.close()

    # yaml.safe_load loads numerical zero values as "None". below is a fix
    if "min_expr" in config_params and config_params["min_expr"] is None:
        config_params["min_expr"] = 0.0
    if "min_prop" in config_params and config_params["min_prop"] is None:
        config_params["min_prop"] = 0.0
    if "padj_thresh" in config_params and config_params["padj_thresh"] is None:
        config_params["padj_thresh"] = 0.0

    return config_params


def check_config():
    '''get config params in a dict then validate the params'''

    err_msg = ""

    config_params = parse_config()
    if isinstance(config_params, str):
        err_msg = config_params
    else:
        err_msg = helpers.validate_parameters(config_params)

    if err_msg:
        delete_user_file("config.yml")

    return err_msg

def get_tsv_rows(filename):
    '''
    returns the rows of the user input file as a 2d array
    '''

    tsv_file = read_user_file(filename)
    data_reader = csv.reader(tsv_file, delimiter="\t")
    rows = list(data_reader)
    tsv_file.close()

    return rows


cleanup_old_sessions()
