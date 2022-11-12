from flask import Blueprint, render_template, request

from app.models.task import Task
from app.blueprints.common import require_valid_task_id

preprocessing_bp = Blueprint('preprocessing_bp', __name__)


@preprocessing_bp.route("/preprocessing")
def preprocessing():
    """Load preprocessing page"""

    task_id = request.args.get("task_id")

    if not task_id or Task.get(task_id) is None:
        task_id = Task.create("preprocessing")

    return render_template("preprocessing.html", title="Preprocessing", task_id=task_id)


@require_valid_task_id
@preprocessing_bp.route("/submit-preprocessing", methods=["POST"])
def submit_preprocessing():
    """Submit preprocessing task"""

    task_id = request.form.get("task_id")
    status_msg = Task.submit(task_id)

    return status_msg


@require_valid_task_id
@preprocessing_bp.route("/confirm-preprocessing-submission", methods=["POST"])
def confirm_preprocessing_submission():
    """Validate submitted datasets and display datasets to load before submission"""

    data_source = request.form.get("source")
    dsets = request.form.get("dsets")
    task_id = request.form.get("task_id")

    config = dict(data_source=data_source, dsets=dsets)

    # Configure task and get confirmation message
    confirmation_message = Task.configure(task_id, config)["status"]

    return confirmation_message
