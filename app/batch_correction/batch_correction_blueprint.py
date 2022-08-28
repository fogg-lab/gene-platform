import os
import time
from flask import Blueprint, render_template, request, session, jsonify, send_from_directory
from ..common import common_blueprint as common
from .. import helpers
from . import bc

# Set the current working directory and relative path to the user files
SCRIPT_PATH = os.path.realpath(__file__)
SCRIPT_DIR = "/".join(SCRIPT_PATH.split("/")[:-1])
USER_FILES_LOCATION = "../user_files"

batch_correction_bp = Blueprint('batch_correction_bp', __name__,
    template_folder = '../templates', static_folder='../static')


@batch_correction_bp.route("/batch-correction")
def batchcorrection():
    """batch correction input form"""

    common.ensure_session_dir()

    cur_uploads, all_uploads = common.list_user_files()

    return render_template("batchcorrection.html", \
        cur_uploads=cur_uploads, all_uploads=all_uploads, title="Batch Correction")


@batch_correction_bp.route("/batch-upload", methods=["POST"])
def batchupload():
    """
    handles uploading counts and coldata files for batch correction
    one file per request
    file contents are in request.data (a bytes object)
    """

    result = {}

    user_filename = request.args.get("user_filename")
    standard_filename = request.headers.get('X_FILENAME')

    if standard_filename not in ["counts.tsv", "coldata.tsv", "filter.txt", \
        "config.yml"]:
        result["error"] = "Unrecognized file."
        return jsonify(result)

    common.save_temp_file(request.data, standard_filename, user_filename)

    if standard_filename == "coldata.tsv":
        result["error_status"] = check_bc_coldata()

    return jsonify(result)


@batch_correction_bp.route("/submit-batch-correction", methods=["POST"])
def submit_batch_correction():
    """Submit a batch correction job"""

    user_dir = common.get_session_dir()

    datatype = request.form.get("data_type")
    reference_level = request.form.get("reference_level")
    contrast_level = request.form.get("contrast_level")
    userdir = session["user_session_dir"]

    expected_bc_counts_path = os.path.join(user_dir, "counts_bc.tsv")

    # Delete any previous batch correction results
    if os.path.isfile(expected_bc_counts_path):
        os.remove(expected_bc_counts_path)

    status_msg = bc.call_bc(userdir, datatype, reference_level, contrast_level)
    if not status_msg:
        status_msg = "Batch correction complete."

    is_output = False
    while not is_output:
        is_output = os.path.exists(expected_bc_counts_path)
        if not is_output:
            time.sleep(0.25)

    return status_msg


@batch_correction_bp.route("/get-batch-corrected-counts")
def get_batch_correction_counts():
    """Returns batch correction results to the client"""

    rel_user_dir = common.get_session_dir()
    abs_user_dir = os.path.abspath(rel_user_dir)
    return send_from_directory(abs_user_dir, "counts_bc.tsv")


def check_bc_coldata():
    """Ensures coldata has batches"""

    coldata = common.get_tsv_rows("coldata.tsv")
    err_msg = helpers.ensure_batches(coldata)

    if err_msg:
        helpers.delete_user_file("coldata.tsv",  common.get_session_dir())

    return err_msg
