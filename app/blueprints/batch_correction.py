import os
import time
from flask import (Blueprint, render_template, request, session, jsonify,
                   send_from_directory)
from app.models.job import Job
from app.job_runner import prepare_job

batch_correction_bp = Blueprint("batch_correction_bp", __name__)


@batch_correction_bp.route("/batch-correction")
def batchcorrection():
    """batch correction input form"""

    common.ensure_session_dir()

    cur_uploads, all_uploads = common.list_user_files()

    return render_template("batchcorrection.html", cur_uploads=cur_uploads,
                           all_uploads=all_uploads, title="Batch Correction")


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

    if standard_filename not in ["counts.tsv", "coldata.tsv", "filter.txt", "config.yml"]:
        result["error"] = "Unrecognized file."
        return jsonify(result)

    common.save_temp_file(request.data, standard_filename, user_filename)

    if standard_filename == "coldata.tsv":
        result["error_status"] = check_batch_correction_coldata()

    return jsonify(result)


@batch_correction_bp.route("/submit-batch-correction", methods=["POST"])
def submit_batch_correction():
    """Submit a batch correction job"""

    user_dir = common.get_session_dir()

    datatype = request.form.get("data_type")
    reference_level = request.form.get("reference_level")
    contrast_level = request.form.get("contrast_level")
    userdir = session["user_session_dir"]

    # Validate input files before calling batch correction
    error_msg = bc_prep.check_bc_input_files(user_dir)
    if error_msg:
        return error_msg

    # Ensure that counts and coldata have the same samples in the same order
    counts_colnames = common.get_tsv_rows("counts.tsv")[0]
    error_msg = helpers.check_coldata_matches_counts(
        counts_colnames, common.get_tsv_rows("coldata.tsv"))
    if error_msg:
        return error_msg

    # Ensure factor levels are present in coldata
    error_msg = helpers.check_factor_levels(reference_level, contrast_level,
                                            common.get_tsv_rows("coldata.tsv"))
    if error_msg:
        return error_msg

    ### Start batch correction job ###
    expected_bc_counts_path = os.path.join(user_dir, "counts_bc.tsv")

    # Delete any previous batch correction results
    if os.path.isfile(expected_bc_counts_path):
        os.remove(expected_bc_counts_path)

    bc_prep.call_batch_correction(userdir, datatype, reference_level, contrast_level)

    is_output = False
    while not is_output:
        is_output = os.path.exists(expected_bc_counts_path)
        if not is_output:
            time.sleep(0.25)

    return "Batch correction complete."


@batch_correction_bp.route("/get-batch-corrected-counts")
def get_batch_correction_counts():
    """Returns batch correction results to the client"""

    rel_user_dir = common.get_session_dir()
    abs_user_dir = os.path.abspath(rel_user_dir)
    return send_from_directory(abs_user_dir, "counts_bc.tsv")
