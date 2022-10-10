from flask import Blueprint, render_template, request, Response
from app.models.job import Job

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

    Job.delete_input_file(job_id, filename)

    return f"{filename} upload cancelled"


@common_bp.route("/get-progress")
def get_console_output():
    """Returns status and updated log of a running job."""

    job_id = request.args.get("job_id")
    last_log_update_line_number = request.args.get("last_log_update_line_number")

    job = Job.get(job_id)

    if job is None:
        return Response("Job not found", status=404)

    log_content = job.get_log_update(last_log_update_line_number)

    return Response(log_content, mimetype='text/plain')
