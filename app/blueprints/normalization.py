from pathlib import Path
from flask import Blueprint, render_template, request, send_from_directory

from app.models.task import Task
from app.blueprints.common import require_valid_task_id

normalization_bp = Blueprint('normalization_bp', __name__)


@normalization_bp.route("/normalization")
def normalization():
    """normalization input form"""

    task_id = request.args.get("task_id")

    if not task_id:
        task_id = Task.create("normalization")

    uploads = Task.list_input_files(task_id)

    return render_template("normalization.html", uploaded_input_files=uploads,
                           title="Normalization")


@require_valid_task_id
@normalization_bp.route("/submit-normalization", methods=["POST"])
def submit_normalization():
    """Submit a normalization task"""

    method = request.form.get("method")
    task_id = request.form.get("task_id")

    Task.configure(task_id, {"method": method})

    status_msg = Task.submit(task_id)

    return status_msg


@require_valid_task_id
@normalization_bp.route("/get-normalized-counts")
def get_normalized_counts():
    """Return the normalized counts file to the client"""

    task_id = request.args.get("task_id")

    filename = "counts_normalized.tsv"

    filepath = Task.get_output_filepath(task_id, "")
    if filepath == "":
        return f"'{filename}' not found.", 204

    filedir = Path(filepath).parent

    return send_from_directory(filedir, filename)
