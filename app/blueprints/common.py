import os
from flask import Blueprint, render_template, request, Response
from app.models.job import Job
from app.job_runner.job_runner import get_job_log_update

common_bp = Blueprint('common_bp', __name__)


@common_bp.route("/")
def index():
    """Load the homepage."""
    return render_template("home.html", title="Welcome!")


@common_bp.route("/cancel-upload", methods=["POST"])
def cancelupload():
    """Remove uploaded file in job directory"""

    filename = request.form.get("filename")
    job_id = request.form.get("job_id")

    job_dir = Job.get_dir(job_id)
    file_path = os.path.join(job_dir, filename)
    if os.path.exists(file_path):
        os.remove(file_path)

    return f"{filename} upload cancelled"


@common_bp.route("/get-console-output")
def get_console_output():
    """Returns the contents of the log file for a job."""

    job_id = request.args.get("job_id")

    job_dir = Job.get_dir(job_id)
    log_path = os.path.join(job_dir, ".log")

    if not os.path.isfile(log_path):
        return ("", 204)

    log_content = get_job_log_update(job_dir)

    return Response(log_content, mimetype='text/plain')
