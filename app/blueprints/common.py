from functools import wraps
from flask import Blueprint, render_template, request, Response

from app.models.task import Task

common_bp = Blueprint('common_bp', __name__)


@common_bp.route("/")
def index():
    """Load the homepage."""
    return render_template("home.html", title="Welcome!")


@common_bp.route("/cancel-upload", methods=["POST"])
def cancelupload():
    """Remove uploaded file in task directory"""

    filename = request.form.get("filename")
    task_id = request.form.get("task_id")

    Task.delete_input_file(task_id, filename)

    return f"{filename} upload cancelled"


@common_bp.route("/get-progress")
def get_console_output():
    """Returns status and updated log of a running task."""

    task_id = request.args.get("task_id")
    last_log_update_line_number = request.args.get("last_log_update_line_number")

    task = Task.get(task_id)

    if task is None:
        return Response("Task not found", status=404)

    log_content = task.get_log_update(last_log_update_line_number)

    return Response(log_content, mimetype='text/plain')


def require_valid_task_id(task_route):
    """Decorator to check if task_id is present in request args and is valid"""
    @wraps(task_route)

    def check_task_id(*args, **kwargs):
        task_id = request.args.get("task_id")

        if Task.get(task_id) is None:
            err_msg = f"Invalid task ID: {task_id}" if task_id else "No task ID provided"
            return 404, err_msg

        return task_route(*args, **kwargs)

    return check_task_id
