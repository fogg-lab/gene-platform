import os
import subprocess
import shutil
import csv
import tempfile
from flask import Blueprint, render_template, request, session, Response, current_app
from flask_login import current_user

jobs_bp = Blueprint('jobs_bp', __name__, template_folder='templates')


@jobs_bp.route("/cancel-upload", methods=["POST"])
def cancelupload():
    """removes uploaded file from user session directory"""

    filename = request.form.get("filename")
    helpers.delete_user_file(filename, get_session_dir())

    return f"{filename} upload cancelled"


@jobs_bp.route("/get-console-output")
def get_console_output():
    """
    returns the contents of the log file in user session directory
    the log file contains terminal output from the analysis script
    """

    log_content = ""
    log_path = os.path.join(get_session_dir(), ".log")

    if os.path.isfile(log_path):
        with open(log_path, "r", encoding="utf-8") as log_file:
            log_content = log_file.read()[-1000:]
        with open(log_path, 'r+', encoding="utf-8") as log_file:
            log_file.truncate(0)
    else:
        return ("", 204)

    return Response(log_content, mimetype='text/plain')
