import os
import time
import base64
import fitz
from flask import Blueprint, render_template, request, session, jsonify
from app.blueprints.correlation.utils import correlation_prep
from app.blueprints.common import common
import app.helpers as helpers

correlation_bp = Blueprint('correlation_bp', __name__)


@correlation_bp.route("/rnaseq-correlation")
def rnaseq_correlation():
    """RNAseq sample correlation page"""

    common.ensure_session_dir()

    cur_uploads, all_uploads = common.list_user_files()

    return render_template("rnaseq_correlation.html",
                           cur_uploads=cur_uploads,
                           all_uploads=all_uploads,
                           title="RNAseq Sample Correlation")


@correlation_bp.route("/upload-rnaseq-correlation", methods=["POST"])
def upload_rnaseq_correlation():
    """handles uploading counts for rnaseq sample correlation"""

    result = {}

    user_filename = request.args.get("user_filename")
    common.save_temp_file(request.data, "counts.tsv", user_filename)

    return jsonify(result)


@correlation_bp.route("/get-correlation-plot")
def get_correlation_plot():
    """Returns specified correlation plot to client"""

    user_dir = common.get_session_dir()
    print(request.args)
    print(request.form)
    corr_method = request.args.get("corr_method")
    img_path = os.path.join(user_dir, f"{corr_method}.png")

    if not os.path.isfile(img_path):
        print(f"{img_path} not found")
        return ('', 204)

    with open(f'{img_path}', 'rb') as img_fp:
        img_data = base64.b64encode(img_fp.read()).decode("utf-8")
        return img_data


@correlation_bp.route("/submit_rnaseq_correlation", methods=["POST"])
def submit_rnaseq_correlation():
    """Submit job to get rnaseq sample correlation plots"""

    user_dir = common.get_session_dir()

    corr_method = request.form.get("corr_method")

    expect_spearman = corr_method != "pearson"
    expect_pearson = corr_method != "spearman"

    expected_pearson_path = os.path.join(user_dir, "pearson.pdf")
    expected_spearman_path = os.path.join(user_dir, "spearman.pdf")

    # Delete any previous correlation results
    helpers.delete_user_file("pearson.pdf", common.get_session_dir())
    helpers.delete_user_file("spearman.pdf", common.get_session_dir())
    helpers.delete_user_file("pearson.png", common.get_session_dir())
    helpers.delete_user_file("spearman.png", common.get_session_dir())

    status_msg = correlation_prep.call_corr(user_dir, corr_method)
    if not status_msg:
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

    helpers.delete_user_file("pearson.pdf", common.get_session_dir())
    helpers.delete_user_file("spearman.pdf", common.get_session_dir())

    return status_msg
