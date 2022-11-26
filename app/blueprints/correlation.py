import os
import base64
from flask import Blueprint, render_template, request

from app.models.task import Task
from app.blueprints.common import require_valid_task_id

correlation_bp = Blueprint('correlation_bp', __name__)


@correlation_bp.route("/rnaseq-correlation")
def rnaseq_correlation():
    """Loads the RNAseq sample correlation page."""

    task_id = request.args.get("task_id")
    if not task_id or Task.get(task_id) is None:
        task_id = Task.create("correlation")

    uploads = Task.list_input_files(task_id)

    return render_template("rnaseq_correlation.html", uploaded_input_files=uploads,
                           task_id=task_id, title="RNAseq Sample Correlation")


@correlation_bp.route("/get-correlation-plot")
@require_valid_task_id
def get_correlation_plot():
    """Returns specified correlation plot to client."""

    task_id = request.args.get("task_id")
    corr_method = request.args.get("corr_method")

    image_filename = f"{corr_method}.png"
    img_path = Task.get_output_filepath(task_id, image_filename)

    if not os.path.isfile(img_path):
        return f"{image_filename} was not found", 204

    with open(f'{img_path}', 'rb') as img_fp:
        img_data = base64.b64encode(img_fp.read()).decode("utf-8")
        return img_data


@correlation_bp.route("/submit-rnaseq-correlation", methods=["POST"])
@require_valid_task_id
def submit_rnaseq_correlation():
    """Submit task to get rnaseq sample correlation plots"""

    corr_method = request.form.get("corr_method")
    task_id = request.form.get("task_id")

    Task.configure(task_id, {"corr_method": corr_method})

    status_msg = Task.submit(task_id)

    return status_msg
