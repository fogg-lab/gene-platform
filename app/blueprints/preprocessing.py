from flask import Blueprint, render_template, request, current_app
from pathlib import Path
import json

from app.models.task import Task
from app.blueprints.common import require_valid_task_id

preprocessing_bp = Blueprint('preprocessing_bp', __name__)


@preprocessing_bp.route("/preprocessing")
def preprocessing():
    """Load preprocessing page"""

    task_id = request.args.get("task_id")

    if not task_id or Task.get(task_id) is None:
        task_id = Task.create("preprocessing")

    gdc_projects_json = Path(current_app.config["GDC_DATA_PATH"]) / "gdc_projects.json"
    gdc_projects_list = []
    with open(gdc_projects_json, "r", encoding="utf-8") as f:
        projects = json.load(f)

    for project_base in [*projects.keys()]:
        for ext in projects[project_base]:
            gdc_projects_list.append(f"{project_base}-{ext}")

    return render_template("preprocessing.html", title="Data Retrieval (GEO/GDC)", task_id=task_id,
                           gdc_projects=gdc_projects_list)


@preprocessing_bp.route("/submit-preprocessing", methods=["POST"])
@require_valid_task_id
def submit_preprocessing():
    """Submit preprocessing task"""

    task_id = request.form.get("task_id")
    status_msg = Task.submit(task_id)

    return status_msg


@preprocessing_bp.route("/confirm-preprocessing-submission", methods=["POST"])
@require_valid_task_id
def confirm_preprocessing_submission():
    """Validate submitted datasets and display datasets to load before submission"""

    data_source = request.form.get("source")
    dsets = request.form.get("dsets")
    task_id = request.form.get("task_id")

    config = dict(data_source=data_source, dsets=dsets)

    # Configure task and get confirmation message
    confirmation_message = Task.configure(task_id, config)["status"]

    return confirmation_message
