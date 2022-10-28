from flask import Blueprint, render_template, request, jsonify, send_from_directory

from app.models.task import Task
from app.blueprints.common import require_valid_task_id

normalization_bp = Blueprint('normalization_bp', __name__, template_folder='templates',
                             static_folder='../../static')

NORMALIZATION_SCRIPT = "Rscript ../../../rscripts/normalize.r"


@normalization_bp.route("/normalization")
def normalization():
    """normalization input form"""

    task_id = request.args.get("task_id")

    if not task_id:
        task_id = Task.create("normalization")

    uploads = Task.list_input_files(task_id)

    return render_template("normalization.html", uploaded_input_files=uploads,
                           title="Normalization")


@require_valid_task_id
@normalization_bp.route("/normalization-upload", methods=["POST"])
def normalization_upload():
    """
    handles uploading counts and coldata files for normalization
    one file per request
    file contents are in request.data (a bytes object)
    """

    user_filename = request.args.get("user_filename")
    standard_filename = request.headers.get('X_FILENAME')
    task_id = request.form.get("task_id")

    result = {}

    if standard_filename not in ("counts.tsv", "coldata.tsv"):
        result["error"] = "Unrecognized file."
        return jsonify(result)

    result = Task.add_input_file(task_id, request.data, standard_filename, user_filename)

    return jsonify(result)


@require_valid_task_id
@normalization_bp.route("/submit-normalization", methods=["POST"])
def submit_normalization():
    """Submit a normalization task"""

    method = request.form.get("method")
    task_id = request.form.get("task_id")

    #expected_norm_counts_path = os.path.join(user_dir, "counts_normalized.tsv")

    # Delete any previous normalization results
    #if os.path.isfile(expected_norm_counts_path):
        #os.remove(expected_norm_counts_path)

    #subprocess.Popen([f"{NORMALIZATION_SCRIPT} {user_dir} {method}"], shell=True)
    #status_msg = "Normalization complete."

    #is_output = False
    #while not is_output:
        #is_output = os.path.exists(expected_norm_counts_path)
        #if not is_output:
            #time.sleep(0.25)

    #return status_msg
    return ""


@require_valid_task_id
@normalization_bp.route("/get-normalized-counts")
def get_normalized_counts():
    """Return the normalized counts file to the client"""

    #rel_user_dir = common.Task.get_dir(task_id)
    #abs_user_dir = os.path.abspath(rel_user_dir)
    #return send_from_directory(abs_user_dir, "counts_normalized.tsv")
    return 404, "Not implemented"
