from pathlib import Path
from flask import Blueprint, render_template, request, send_from_directory, redirect, url_for
from app.models.job import Job
from app.exceptions import InvalidJobOutput

preprocessing_bp = Blueprint('preprocessing_bp', __name__)

PREP_GDC_SCRIPT = "prep_gdc.r"
PREP_GEO_SCRIPT = "prep_geo.r"
JOB_TYPE = "preprocessing"


@preprocessing_bp.route("/preprocessing")
def preprocessing():
    """Load preprocessing page"""

    job_id = request.args.get("job_id")

    job = Job.get(job_id)
    if job is None:
        job = Job.create(JOB_TYPE)

    uploads = job.list_input_files()

    return render_template(
        "preprocessing.html",
        all_uploads=uploads, title="Preprocessing"
    )


@preprocessing_bp.route("/submit-preprocessing", methods=["POST"])
def submit_preprocessing():
    """Submit preprocessing job"""

    job_id = request.args.get("job_id")
    job = Job.get(job_id)

    if job is None:
        return redirect(url_for("preprocessing_bp.preprocessing"))

    status_msg = job.submit()

    return status_msg


@preprocessing_bp.route("/confirm-preprocessing-submission", methods=["POST"])
def confirm_preprocessing_submission():
    """Validate submitted datasets and display datasets to load before submission"""

    data_source = request.form.get("source")
    dsets = request.form.get("dsets")
    job_id = request.args.get("job_id")

    job = Job.get(job_id)

    if job is None:
        return redirect(url_for("preprocessing_bp.preprocessing"))

    config = dict(data_source=data_source, dsets=dsets)

    # Configure job and get confirmation message
    confirmation_message = job.configure(config)

    return confirmation_message


@preprocessing_bp.route("/get-preprocessed-data")
def get_preprocessed_data():
    """Return preprocessed data for user to download"""

    job_id = request.args.get("job_id")

    job = Job.get(job_id)

    if job is None:
        return redirect(url_for("preprocessing_bp.preprocessing"))

    processed_data = Job.get_output()
    if len(processed_data) != 1 or not processed_data[0].endswith(".zip"):
        raise InvalidJobOutput("The final output for a preprocessing job must be "
                               "a single zip file.")

    processed_zip_path = Path(processed_data[0])
    output_dir = processed_zip_path.parent
    output_fname = processed_zip_path.name

    return send_from_directory(output_dir, output_fname)
