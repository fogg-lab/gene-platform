import os
from flask import (Blueprint, render_template, request, jsonify,
                   send_from_directory)
from app.models.job import Job
from app.job_utils.job_runner import add_input_file, list_input_files

batch_correction_bp = Blueprint("batch_correction_bp", __name__)


@batch_correction_bp.route("/batch-correction")
def batchcorrection():
    """batch correction input form"""

    job_id = request.form.get("job_id")
    job_dir = Job.get_dir(job_id)

    cur_uploads, all_uploads = list_input_files(job_dir)

    return render_template("batchcorrection.html", cur_uploads=cur_uploads,
                           all_uploads=all_uploads, title="Batch Correction")


@batch_correction_bp.route("/batch-upload", methods=["POST"])
def batchupload():
    """
    handles uploading counts and coldata files for batch correction
    one file per request
    file contents are in request.data (a bytes object)
    """

    job_id = request.form.get("job_id")
    job_dir = Job.get_dir(job_id)

    result = {}

    user_filename = request.args.get("user_filename")
    standard_filename = request.headers.get('X_FILENAME')

    """ OLD CODE:
    if standard_filename not in ["counts.tsv", "coldata.tsv", "filter.txt", "config.yml"]:
        result["error"] = "Unrecognized file."
        return jsonify(result)
    if standard_filename == "coldata.tsv":
        result["error_status"] = check_batch_correction_coldata()
    """

    result = add_input_file(request.data, job_dir, standard_filename,
                            user_filename, "batch_correction")
    return jsonify(result)


@batch_correction_bp.route("/submit-batch-correction", methods=["POST"])
def submit_batch_correction():
    """Submit a batch correction job"""

    job_id = request.form.get("job_id")
    job_dir = Job.get_dir(job_id)

    datatype = request.form.get("data_type")
    reference_level = request.form.get("reference_level")
    contrast_level = request.form.get("contrast_level")

    return "Batch correction complete."


@batch_correction_bp.route("/get-batch-corrected-counts")
def get_batch_correction_counts():
    """Returns batch correction results to the client"""

    job_id = request.form.get("job_id")
    job_dir = Job.get_dir(job_id)

    abs_user_dir = os.path.abspath(rel_user_dir)
    return send_from_directory(abs_user_dir, "counts_bc.tsv")
