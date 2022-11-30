import sys
from functools import wraps
from pathlib import Path
from flask import (Blueprint, render_template, request, jsonify, Response,
                   send_from_directory, redirect, url_for)

from app.models.task import Task

common_bp = Blueprint('common_bp', __name__)


@common_bp.route("/")
def index():
    """Load the homepage."""
    return render_template("home.html", title="Welcome!")


def require_valid_task_id(task_route):
    """Decorator to check if task_id is present in request args and is valid"""
    @wraps(task_route)

    def check_task_id(*args, **kwargs):
        task_id = request.form.get("task_id")
        if not task_id:
            task_id = request.args.get("task_id")

        if not task_id or Task.get(task_id) is None:
            err_msg = f"Invalid task ID: {task_id}" if task_id else "No task ID provided"
            return err_msg, 204

        return task_route(*args, **kwargs)

    return check_task_id


@common_bp.route("/upload", methods=["POST"])
@require_valid_task_id
def upload():
    """Receive an uploaded file for a task."""

    user_fname = request.args.get("user_filename")
    standard_fname = request.headers.get('X_FILENAME')
    task_id = request.args.get("task_id")

    print(f"User filename: {user_fname}", file=sys.stderr)

    result = Task.add_input_file(task_id, request.data, standard_fname, user_fname)

    return jsonify(result)


@common_bp.route("/cancel-upload", methods=["POST"])
@require_valid_task_id
def cancelupload():
    """Remove uploaded file in task directory"""

    filename = request.form.get("filename")
    task_id = request.form.get("task_id")

    Task.delete_input_file(task_id, filename)

    return f"{filename} upload cancelled"


@common_bp.route("/get-progress")
@require_valid_task_id
def get_console_output():
    """Returns status and updated log of a running task."""

    last_log_offset = int(request.args.get("last_log_offset"))

    task_id = request.args.get("task_id")

    log_content, last_log_offset = Task.get_log_update(task_id, last_log_offset)
    task_status = Task.get_status(task_id)

    # Append metadata to the log content (last log offset and task status)
    log_content += ('\n###PROGRESS_METADATA:{'
        f'"last_log_offset":{last_log_offset},'
        f'"task_status":"{task_status}"'
        '}')

    return Response(log_content, mimetype='text/plain')


@common_bp.route("/tasks")
def tasks():
    """Load the tasks page populated with all tasks associated with the current user."""

    result = []
    all_tasks = Task.get_user_tasks()
    for task in all_tasks:
        task_data = {}
        task_data["task_id"] = task.task_id
        task_data['task_type'] = task.task_type
        task_data['status'] = task.status
        task_data['created_at'] = task.created_at
        task_data['updated_at'] = task.updated_at

        result.append(task_data)

    result = sorted(result, key=lambda d: d['updated_at'], reverse=True)

    return render_template("tasks.html", title="Tasks", user_tasks=result)


@common_bp.route("/download-task-output/<task_id>")
def get_task_output(task_id):
    """Return a compressed zip file of the task output to the client."""

    if not task_id or Task.get(task_id) is None:
        return f"Invalid task ID: {task_id}" if task_id else "No task ID provided", 204

    zip_out_path = Task.create_task_zip(task_id)

    if not zip_out_path:
        return "No output files to download", 204

    out_dir = Path(zip_out_path).parent
    out_filename = Path(zip_out_path).name

    return send_from_directory(out_dir, out_filename, as_attachment=True)


@common_bp.route("/remove-task/<task_id>")
def remove_task(task_id):
    """Remove a task from the database and delete its directory."""

    if not task_id or Task.get(task_id) is None:
        return f"Invalid task ID: {task_id}" if task_id else "No task ID provided", 204

    Task.delete(task_id)

    return redirect(url_for('common_bp.tasks'))
