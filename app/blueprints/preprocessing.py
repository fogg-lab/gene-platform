from pathlib import Path
from flask import Blueprint, render_template, request, send_from_directory

from app.models.task import Task
from app.blueprints.common import require_valid_task_id

preprocessing_bp = Blueprint('preprocessing_bp', __name__)


@preprocessing_bp.route("/preprocessing")
def preprocessing():
    """Load preprocessing page"""

    task_id = request.args.get("task_id")

    if Task.get(task_id) is None:
        task_id = Task.create("preprocessing")

    return render_template("preprocessing.html", title="Preprocessing")


@require_valid_task_id
@preprocessing_bp.route("/submit-preprocessing", methods=["POST"])
def submit_preprocessing():
    """Submit preprocessing task"""

    task_id = request.args.get("task_id")
    status_msg = Task.submit(task_id)

    return status_msg


@require_valid_task_id
@preprocessing_bp.route("/confirm-preprocessing-submission", methods=["POST"])
def confirm_preprocessing_submission():
    """Validate submitted datasets and display datasets to load before submission"""

    data_source = request.form.get("source")
    dsets = request.form.get("dsets")
    task_id = request.args.get("task_id")

    config = dict(data_source=data_source, dsets=dsets)

    # Configure task and get confirmation message
    confirmation_message = Task.configure(task_id, config)

    return confirmation_message


@require_valid_task_id
@preprocessing_bp.route("/get-preprocessed-data")
def get_preprocessed_data():
    """Return preprocessed data for user to download"""

    task_id = request.args.get("task_id")

    processed_data = Task.create_task_zip(task_id)
    processed_zip_path = Path(processed_data)

    output_dir = processed_zip_path.parent
    output_fname = processed_zip_path.name

    return send_from_directory(output_dir, output_fname)
