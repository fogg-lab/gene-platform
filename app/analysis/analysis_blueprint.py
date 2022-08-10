import os
import subprocess
import csv
import yaml
import time
from ..common import common_blueprint as common
from .. import helpers
from flask import Blueprint, render_template, request, redirect, url_for, \
    session, Response, jsonify, send_from_directory

RNA_SEQ_SCRIPT = "Rscript ../../rscripts/rnaseq.R"
MICROARRAY_SCRIPT = "Rscript ../../rscripts/microarray.R"

# Set the current working directory and relative path to the user files
SCRIPT_PATH = os.path.realpath(__file__)
SCRIPT_DIR = "/".join(SCRIPT_PATH.split("/")[:-1])
USER_FILES_LOCATION = "../user_files"

os.chdir(SCRIPT_DIR)

analysis_bp = Blueprint('analysis_bp', __name__,
    template_folder='../templates', static_folder='../static')


@analysis_bp.route("/setup")
def setup():
    '''first input page (file uploads)'''

    common.ensure_session_dir()

    cur_uploads, all_uploads = common.list_user_files()

    return render_template("uploads_form.html", \
        cur_uploads=cur_uploads, all_uploads=all_uploads, title="Uploads")


@analysis_bp.route("/upload", methods=["POST"])
def upload():
    '''
    handles uploading counts, coldata, filter, and config file
    one file per request
    file contents are in request.data (a bytes object)
    '''

    result = {}

    user_filename = request.args.get("user_filename")
    standard_filename = request.headers.get('X_FILENAME')

    if standard_filename not in ["counts.tsv", "coldata.tsv", "filter.txt", \
        "config.yml"]:
        result["error"] = "Unrecognized file."
        return jsonify(result)

    common.save_temp_file(request.data, standard_filename, user_filename)

    if standard_filename == "config.yml":
        result["error_status"] = check_config()

    return jsonify(result)


@analysis_bp.route("/parameters", methods=["GET"])
def parameters():
    '''loads the parameter form'''

    params = parse_config()
    return render_template("parameters_form.html", params=params,\
         title="Parameters")


@analysis_bp.route("/submit", methods=["POST"])
def submit():
    '''
    submit_parameters()
    submits input for analysis if valid
    generates a config file from the parameter form
    '''

    # get whether analysis is microarray or RNA-Seq
    data_type = request.form.get("data_type")

    # remove old output
    helpers.delete_user_file("filter_output.tsv", common.get_session_dir())
    helpers.delete_user_file("output.tsv", common.get_session_dir())

    for filename in os.listdir(session["user_session_dir"]):
        if filename.endswith(".png"):
            helpers.delete_user_file(filename, common.get_session_dir())

    call_analysis(data_type)

    wait_for_output()

    return "analysis completed"


@analysis_bp.route("/confirm-analysis-submission", methods=["POST"])
def confirm_analysis_submission():
    '''validate input and display formula before submission'''

    # get whether analysis is microarray or RNA-Seq, and get params
    data_type = request.form.get("data_type")
    params = helpers.get_request_parameters(request.form, data_type)

    # generate config file from the form parameters
    generate_config(params)
    
    # validate the config, counts and coldata
    config_file_error = check_config()

    counts_colnames = common.get_tsv_rows("counts.tsv")[0]

    coldata_counts_match_error = helpers.check_coldata_matches_counts(
        counts_colnames, common.get_tsv_rows("coldata.tsv"))
    
    factor_levels_error = helpers.check_factor_levels(
        params, common.get_tsv_rows("coldata.tsv"))

    # get the analysis formula to display for the user
    confirmation_message = ""

    if coldata_counts_match_error:
        confirmation_message += f"<p>Error: {coldata_counts_match_error}</p>"
    if factor_levels_error:
        confirmation_message += f"<p>Error: {factor_levels_error}</p>"
    if config_file_error:
        confirmation_message += f"<p>Error: {config_file_error}</p>"

    if not confirmation_message:
        confirmation_message = helpers.get_analysis_confirmation_msg(params)

    return confirmation_message


@analysis_bp.route("/display")
def display_output():
    '''
    display output in a table
    '''

    # check here if output.tsv exists and errors.txt doesn't
    unfiltered_output_path = f"{session['user_session_dir']}output.tsv"
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


@analysis_bp.route("/plots")
def plots():
    '''
    display plots
    '''

    common.ensure_session_dir()

    session_dir = f"{session['user_session_dir']}"
    plot_filenames = {}

    if session_dir:
        for filename in os.listdir(session_dir):
            if "mean" in filename and ".png" and "unfiltered" in filename:
                plot_filenames["unfiltered_mean_variance"] = filename
            elif "mean" in filename and ".png" in filename:
                plot_filenames["filtered_mean_variance"] = filename
            elif "unfilt" in filename and ".png" in filename:
                plot_filenames["unfiltered_volcano"] = filename
            elif ".png" in filename:
                plot_filenames["filtered_volcano"] = filename

    return render_template("plots.html", plot_filenames=plot_filenames)


@analysis_bp.route('/getplot/<filename>')
def getplot(filename):
    rel_user_dir = common.get_session_dir()
    abs_user_dir = os.path.abspath(rel_user_dir)
    return send_from_directory(abs_user_dir, filename)


@analysis_bp.route("/reset", methods=["GET"])
def reset():
    '''Deletes files from a users session'''

    helpers.delete_user_file("counts.tsv",  common.get_session_dir())
    helpers.delete_user_file("coldata.tsv",  common.get_session_dir())
    helpers.delete_user_file("filter.txt", common.get_session_dir())
    helpers.delete_user_file("config.yml", common.get_session_dir())
    helpers.delete_user_file("filter_output.tsv", common.get_session_dir())
    helpers.delete_user_file("output.tsv", common.get_session_dir())

    return redirect(url_for("common_bp.index"))


@analysis_bp.route("/get-unfiltered-tsv")
def get_unfiltered_tsv():
    '''download unfiltered output'''

    unfiltered_output = open(f"{session['user_session_dir']}output.tsv", \
        encoding="UTF-8")

    return Response(
        unfiltered_output,
        mimetype='text/csv',
        headers={'Content-disposition':
                'attachment; filename=output.tsv'})


@analysis_bp.route("/get-filtered-tsv")
def get_filtered_tsv():
    '''download filtered output'''

    filtered_output = open(f"{session['user_session_dir']}filter_output.tsv",\
         encoding="UTF-8")

    return Response(
        filtered_output,
        mimetype="text/csv",
        headers={"Content-disposition":
                "attachment; filename=filter_output.tsv"})


def call_analysis(data_type):
    '''
    calls microarray or rna-seq analysis depending on data_type
    sends the session_id as an argument
    '''

    log = common.get_session_dir() + "log"

    if data_type == 'microarray':
        subprocess.Popen([f"{MICROARRAY_SCRIPT} {session['session_id']} "
                          f"1> {log} 2>& 1"], shell=True)
    elif data_type == 'RNA-Seq':
        subprocess.Popen([f"{RNA_SEQ_SCRIPT} {session['session_id']} "
                          f"1> {log} 2>& 1"], shell=True)


def wait_for_output():
    '''waits until output shows up in users session directory'''

    unfilt_output_path = f"{session['user_session_dir']}output.tsv"
    filt_output_path = f"{session['user_session_dir']}filter_output.tsv"
    filter_path = f"{session['user_session_dir']}filter.txt"

    is_output = False
    while not is_output:
        is_output = os.path.exists(unfilt_output_path)
        if is_output and os.path.exists(filter_path):
            is_output = os.path.exists(filt_output_path)

        if not is_output:
            time.sleep(0.5)


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


def parse_config():
    ''' parses config into a dict of parameters '''

    session_dir = common.get_session_dir()
    config_file_path = session_dir + "config.yml"
    config_params = {}

    if session_dir and os.path.exists(config_file_path):
        config_file = open(config_file_path, encoding="UTF-8")
        config_params = yaml.load(config_file, Loader=yaml.FullLoader)
        config_file.close()

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
        helpers.delete_user_file("config.yml", common.get_session_dir())

    return err_msg
