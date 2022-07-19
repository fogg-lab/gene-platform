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
PREP_GDC_SCRIPT = "Rscript ../../rscripts/prep_gdc.r"
PREP_GEO_SCRIPT = "Rscript ../../rscripts/prep_geo.r"

os.chdir(SCRIPT_DIR)

preprocessing_bp = Blueprint('preprocessing_bp', __name__, template_folder='../templates', \
    static_folder='../static')

@preprocessing_bp.route("/preprocessing")
def preprocessing():
    '''todo'''

    common.ensure_session_dir()

    cur_uploads, all_uploads = common.list_user_files()

    return render_template("preprocessing.html", \
        cur_uploads=cur_uploads, all_uploads=all_uploads, title="Batch Correction")

@preprocessing_bp.route("/submit_preprocessing", methods=["POST"])
def submit_preprocessing():
    data_source = request.form.get("data_source")
    target = request.form.get("target")
    user_dir = session["user_session_dir"]

    expected_counts_path = f"{user_dir}counts_preprocessed.tsv"
    expected_coldata_path = f"{user_dir}coldata_preprocessed.tsv"

    # Delete any previous preprocessing results
    if os.path.isfile(expected_counts_path):
        os.remove(expected_counts_path)
    if os.path.isfile(expected_coldata_path):
        os.remove(expected_coldata_path)

    if data_source == "gdc":
        subprocess.Popen([f"{PREP_GDC_SCRIPT} {user_dir} {target}"], shell=True)
    elif data_source == "geo":
        subprocess.Popen([f"{PREP_GEO_SCRIPT} {user_dir} {target}"], shell=True)
    
    status_msg = "Preprocessing complete."

    is_output = False
    while not is_output:
        is_output = os.path.exists(expected_counts_path)\
            and os.path.exists(expected_coldata_path)
        if not is_output:
            time.sleep(0.25)

    return status_msg


@preprocessing_bp.route("/get_preprocessed_counts")
def get_preprocessed_counts():
    rel_user_dir = common.get_session_dir()
    abs_user_dir = os.path.abspath(rel_user_dir)
    return send_from_directory(abs_user_dir, "counts_preprocessed.tsv")


@preprocessing_bp.route("/get_preprocessed_coldata")
def get_preprocessed_coldata():
    rel_user_dir = common.get_session_dir()
    abs_user_dir = os.path.abspath(rel_user_dir)
    return send_from_directory(abs_user_dir, "coldata_preprocessed.tsv")
