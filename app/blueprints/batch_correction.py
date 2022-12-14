from pathlib import Path
from flask import Blueprint, render_template, request, send_from_directory

from app.models.task import Task
from app.helper import format_msg_list_html
from app.blueprints.common import require_valid_task_id

batch_correction_bp = Blueprint("batch_correction_bp", __name__)


@batch_correction_bp.route("/batch-correction")
def batchcorrection():
    """Loads the batch correction page."""

    task_id = request.args.get("task_id")

    if not task_id or Task.get(task_id) is None:
        task_id = Task.create("batch correction")

    uploads = Task.list_input_files(task_id)

    return render_template("batchcorrection.html", uploaded_input_files=uploads,
                           task_id=task_id, title="Batch Correction")


@batch_correction_bp.route("/submit-batch-correction", methods=["POST"])
@require_valid_task_id
def submit_batch_correction():
    """Submit a batch correction task"""

    task_id = request.form.get("task_id")

    config = {}

    config["data_type"] = request.form.get("data_type")
    config["reference_level"] = request.form.get("reference_level")
    config["contrast_level"] = request.form.get("contrast_level")

    config_status = Task.configure(task_id, config)
    config_error_msg = format_msg_list_html(config_status.get("errors"), tag="p")

    if config_error_msg:
        return config_error_msg

    task_status = Task.submit(task_id)

    return task_status


@batch_correction_bp.route("/get-batch-corrected-counts")
@require_valid_task_id
def get_batch_correction_counts():
    """Returns batch correction results to the client"""

    task_id = request.args.get("task_id")

    bc_counts_path = Task.get_output_filepath(task_id, "batch_corrected_counts.tsv")
    output_file_dir = Path(bc_counts_path).parent if bc_counts_path else None

    if not output_file_dir:
        return "Batch correction output was not found.", 200

    return send_from_directory(output_file_dir, "batch_corrected_counts.tsv")
