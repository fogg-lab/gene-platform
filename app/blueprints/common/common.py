import os
import subprocess
import shutil
import csv
import tempfile
from flask import Blueprint, render_template, request, session, Response, current_app
import app.helpers as helpers

common_bp = Blueprint('common_bp', __name__, template_folder='templates')


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
    log = read_user_file(".log", get_session_dir())
    if log:
        log_content = log.read()[-1000:]
        log.close()
        clear_user_file(".log", get_session_dir())
    else:
        return (log_content, 204)

    return Response(log_content, mimetype='text/plain')


def save_temp_file(file_contents, standard_filename, user_filename):
    """
    Saves uploaded file to user session directory

    Parameters:
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


def read_user_file(filename, user_dir):
    """
    Opens a user file for reading. Caller is responsible for closing the file
    Parameters:
        filename: string
        user_dir: string
    Returns:
        file object
    """

    user_file = None

    if user_dir:
        filepath = os.path.join(user_dir, filename)
        if os.path.isfile(filepath):
            user_file = open(filepath, "r", encoding="UTF-8")

    return user_file


def clear_user_file(filename, user_dir):
    """
    Clears all content from a user file if it exists
    Parameters:
        filename: string
        user_dir: string
    """
    if user_dir:
        filepath = os.path.join(user_dir, filename)
        if os.path.isfile(filepath):
            f = open(filepath, "w", encoding="UTF-8")
            f.truncate()
            f.seek(0)
            f.close()


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


def get_session_dir():
    """
    Returns path to user's session directory.
    If none exists, ensure_session_dir() is called to create one
    """

    if "user_session_dir" not in session:
        ensure_session_dir()

    return session["user_session_dir"]


def ensure_session_dir():
    """
    if there is no directory for the user session, create one now
    directory path is stored in the user session variable 'user_session_dir'
    """

    if ("user_session_dir" not in session or
        not os.path.exists(session["user_session_dir"])):
        user_files_dir = current_app.config["USER_FILES_LOCATION"]
        temp_dir = tempfile.mkdtemp(dir=user_files_dir)
        os.chmod(temp_dir, 0o777) # give everyone rwx permission for the dir

        session["user_session_dir"] = os.path.abspath(temp_dir)
        session["session_id"] = temp_dir.split("/")[-1:]


def get_tsv_rows(filename):
    """returns the rows of the user input file as a 2d array"""

    user_dir = get_session_dir()

    tsv_file = read_user_file(filename, user_dir)
    data_reader = csv.reader(tsv_file, delimiter="\t")
    rows = list(data_reader)
    tsv_file.close()

    return rows


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
