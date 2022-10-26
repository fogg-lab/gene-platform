from pathlib import Path
from flask import Blueprint, render_template, request, send_from_directory, redirect, url_for
from app.models.task import Task
from app.exceptions import InvalidTaskOutput

preprocessing_bp = Blueprint('preprocessing_bp', __name__)

PREP_GDC_SCRIPT = "prep_gdc.r"
PREP_GEO_SCRIPT = "prep_geo.r"
TASK_TYPE = "preprocessing"


@preprocessing_bp.route("/preprocessing")
def preprocessing():
    """Load preprocessing page"""

    task_id = request.args.get("task_id")

    task = Task.get(task_id)
    if task is None:
        task = Task.create(TASK_TYPE)

    uploads = task.list_input_files()

    return render_template(
        "preprocessing.html",
        all_uploads=uploads, title="Preprocessing"
    )


@preprocessing_bp.route("/submit-preprocessing", methods=["POST"])
def submit_preprocessing():
    """Submit preprocessing task"""

    task_id = request.args.get("task_id")
    task = Task.get(task_id)

    if task is None:
        return redirect(url_for("preprocessing_bp.preprocessing"))

    status_msg = task.submit()

    return status_msg


@preprocessing_bp.route("/confirm-preprocessing-submission", methods=["POST"])
def confirm_preprocessing_submission():
    """Validate submitted datasets and display datasets to load before submission"""

    data_source = request.form.get("source")
    dsets = request.form.get("dsets")
    task_id = request.args.get("task_id")

    task = Task.get(task_id)

    if task is None:
        return redirect(url_for("preprocessing_bp.preprocessing"))

    config = dict(data_source=data_source, dsets=dsets)

    # Configure task and get confirmation message
    confirmation_message = task.configure(config)

    return confirmation_message


@preprocessing_bp.route("/get-preprocessed-data")
def get_preprocessed_data():
    """Return preprocessed data for user to download"""

    task_id = request.args.get("task_id")

    task = Task.get(task_id)

    if task is None:
        return redirect(url_for("preprocessing_bp.preprocessing"))

    processed_data = Task.get_output()
    if len(processed_data) != 1 or not processed_data[0].endswith(".zip"):
        raise InvalidTaskOutput("The final output for a preprocessing task must be "
                               "a single zip file.")

    processed_zip_path = Path(processed_data[0])
    output_dir = processed_zip_path.parent
    output_fname = processed_zip_path.name

    return send_from_directory(output_dir, output_fname)
