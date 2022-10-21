import os
import time
import base64
import fitz
from flask import Blueprint, render_template, request, session, jsonify
from app.models.job import Job

correlation_bp = Blueprint('correlation_bp', __name__)

JOB_TYPE = "correlation"

@correlation_bp.route("/rnaseq-correlation")
def rnaseq_correlation():
    """RNAseq sample correlation page"""

    if not job_id:
        job_id = Job.create(JOB_TYPE)

    uploads = Job.list_input_files(job_dir)

    return render_template("rnaseq_correlation.html", cur_uploads=uploads,
                           title="RNAseq Sample Correlation")


@correlation_bp.route("/upload-rnaseq-correlation", methods=["POST"])
def upload_rnaseq_correlation():
    """handles uploading counts for rnaseq sample correlation"""

    job_id = request.form.get("job_id")
    job_dir = Job.get_dir(job_id)
    save_path = os.path.join(job_dir, "counts.tsv")
    data = request.data

    result = dict()
    result["status_msg"] = add_input_file(data, save_path, "correlation")

    return jsonify(result)


@correlation_bp.route("/get-correlation-plot")
def get_correlation_plot():
    """Returns specified correlation plot to client"""

    job_id = request.args.get("job_id")
    job_dir = Job.get_dir(job_id)
    corr_method = request.args.get("corr_method")
    img_path = os.path.join(job_dir, f"{corr_method}.png")

    if not os.path.isfile(img_path):
        print(f"{img_path} not found")
        return ('', 204)

    with open(f'{img_path}', 'rb') as img_fp:
        img_data = base64.b64encode(img_fp.read()).decode("utf-8")
        return img_data


@correlation_bp.route("/submit-rnaseq-correlation", methods=["POST"])
def submit_rnaseq_correlation():
    """Submit job to get rnaseq sample correlation plots"""

    corr_method = request.form.get("corr_method")
    job_id = request.form.get("job_id")
    job_dir = Job.get_dir(job_id)
    job_output_dir = os.path.join(job_dir, "output")

    expect_spearman = corr_method != "pearson"
    expect_pearson = corr_method != "spearman"

    expected_pearson_path = os.path.join(job_output_dir, "pearson.pdf")
    expected_spearman_path = os.path.join(job_output_dir, "spearman.pdf")

    # Delete any previous correlation results
    for fname in job_output_dir:
        if "pearson" in fname or "spearman" in fname:
            file_path = os.path.join(job_output_dir, fname)
            os.remove(file_path)

    correlation_prep.call_corr(user_dir, corr_method)

    status_msg = "Done computing sample correlations."

    is_output = False
    while not is_output:
        if expect_pearson and expect_spearman:
            is_output = (os.path.isfile(expected_pearson_path)
                         and os.path.isfile(expected_spearman_path))
        elif expect_pearson:
            is_output = os.path.isfile(expected_pearson_path)
        elif expect_spearman:
            is_output = os.path.isfile(expected_spearman_path)
        if not is_output:
            time.sleep(0.2)

    # Convert to png for easier display in browser
    save_path = session['user_session_dir']
    for expected_path in [expected_pearson_path, expected_spearman_path]:
        if os.path.isfile(expected_path):
            print(f"Converting {expected_path} to png\n")
            doc = ""
            while not doc:
                try:
                    doc = fitz.open(expected_path)
                except fitz.fitz.EmptyFileError:
                    time.sleep(0.2)
            for page in doc:
                pix = page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5))

                save_path = f"{expected_path[:-4]}.png"
                pix.save(save_path)

    helpers.delete_user_file("pearson.pdf", common.Job.get_dir(job_id))
    helpers.delete_user_file("spearman.pdf", common.Job.get_dir(job_id))

    return status_msg
