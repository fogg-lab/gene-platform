import os
import subprocess
from ..common import common_blueprint as common
import time
from flask import Blueprint, render_template, request, session, jsonify, \
    send_from_directory, session

# Set the current working directory and relative path to the user files
SCRIPT_PATH = os.path.realpath(__file__)
SCRIPT_DIR = "/".join(SCRIPT_PATH.split("/")[:-1])
USER_FILES_LOCATION = "../user_files"
NORMALIZATION_SCRIPT = "Rscript ../../rscripts/normalize.r"

os.chdir(SCRIPT_DIR)

normalization_bp = Blueprint('normalization_bp', __name__, template_folder='../templates', \
    static_folder='../static')


@normalization_bp.route("/normalization")
def normalization():
    '''normalization input form'''

    common.ensure_session_dir()

    cur_uploads, all_uploads = common.list_user_files()

    return render_template("normalization.html", \
        cur_uploads=cur_uploads, all_uploads=all_uploads, title="Normalization")


@normalization_bp.route("/normalization_upload", methods=["POST"])
def normalization_upload():
    '''
    handles uploading counts and coldata files for normalization
    one file per request
    file contents are in request.data (a bytes object)
    '''

    result = {}

    user_filename = request.args.get("user_filename")
    standard_filename = request.headers.get('X_FILENAME')

    if standard_filename not in {"counts.tsv", "coldata.tsv"}:
        result["error"] = "Unrecognized file."
        return jsonify(result)

    common.save_temp_file(request.data, standard_filename, user_filename)

    return jsonify(result)


@normalization_bp.route("/submit-normalization", methods=["POST"])
def submit_normalization():

    method = request.form.get("method")
    user_dir = session["user_session_dir"]

    expected_norm_counts_path = f"{user_dir}counts_normalized.tsv"

    # Delete any previous normalization results
    if os.path.isfile(expected_norm_counts_path):
        os.remove(expected_norm_counts_path)

    subprocess.Popen([f"{NORMALIZATION_SCRIPT} {user_dir} {method}"], shell=True)
    status_msg = "Normalization complete."

    is_output = False
    while not is_output:
        is_output = os.path.exists(expected_norm_counts_path)
        if not is_output:
            time.sleep(0.25)

    return status_msg


@normalization_bp.route("/get-normalized-counts")
def get_normalized_counts():
    rel_user_dir = common.get_session_dir()
    abs_user_dir = os.path.abspath(rel_user_dir)
    return send_from_directory(abs_user_dir, "counts_normalized.tsv")
