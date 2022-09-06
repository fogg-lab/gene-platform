import os
import subprocess
import shutil
import csv
import tempfile
from flask import Blueprint, render_template, request, session, Response, current_app
from app.models.job import Job
from app.job_runner import prepare_job

common_bp = Blueprint('common_bp', __name__)


@common_bp.route("/")
def index():
    """main page"""

    return render_template("home.html", title="Welcome!")


@common_bp.route("/cancel-upload", methods=["POST"])
def cancelupload():
    """removes uploaded file from user session directory"""

    filename = request.form.get("filename")
    helpers.delete_user_file(filename, get_session_dir())

    return f"{filename} upload cancelled"


@common_bp.route("/get-console-output")
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


def save_temp_file(file_contents, standard_filename, user_filename):
    """
    Saves uploaded file to user session directory

    Args:
        file_contents: bytes object
        filename: string
    """

    # create session directory if none exists
    user_dir = get_session_dir()

    user_file_path = os.path.join(user_dir, standard_filename)

    if os.path.exists(user_file_path):
        os.remove(user_file_path)

    # save the user-specified filename
    session[standard_filename] = user_filename

    user_file = open(user_file_path, "w", encoding="UTF-8")
    lines = file_contents.split(b'\n')

    # copy the file line by line, adding a trailing newline if none exists
    for line in lines:
        user_file.write(f"{line.decode('UTF-8')}\n")

    user_file.close()


def cleanup_old_sessions():
    """Clean up other sessions older than 4 hours (14400 seconds)"""

    user_files_dir = current_app.config["USER_FILES_LOCATION"]

    for old_dir in os.listdir():
        old_dir = os.path.join(user_files_dir, old_dir)

        if old_dir.split("/")[-1] != ".gitkeep":
            get_age = f"$(($(date +%s) - $(date +%s -r {old_dir})))"

            # Run bash command and return the stdout output
            age_seconds = subprocess.Popen([f"echo {get_age}"],
                stdout=subprocess.PIPE, shell=True).communicate()[0]

            try:
                age_seconds = int(age_seconds)
                if age_seconds > 14400:
                    shutil.rmtree(old_dir)
            except ValueError:
                pass


def list_user_files():
    all_files = os.listdir(session["user_session_dir"])
    current_uploads = {}
    all_uploads = {}
    standard_input_files = ["counts", "coldata", "filter", "config"]

    for filename in ["counts.tsv", "coldata.tsv", "filter.txt", "config.yml"]:
        current_uploads[filename] = False
    for filename in all_files:
        if filename in current_uploads:
            current_uploads[filename] = True

    for filename in all_files:
        fname_base = ""
        if '-' in filename:
            fname_base = filename.split("-")[0]
        elif '.' in filename:
            fname_base = filename.split(".")[0]
        else:
            continue
        if fname_base in standard_input_files and filename in session:
            all_uploads[filename] = session[filename]

    return current_uploads, all_uploads


cleanup_old_sessions()
