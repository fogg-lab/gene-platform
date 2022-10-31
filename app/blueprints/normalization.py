from flask import Blueprint, render_template, request, send_from_directory

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
@normalization_bp.route("/submit-normalization", methods=["POST"])
def submit_normalization():
    """Submit a normalization task"""

    method = request.form.get("method")
    task_id = request.form.get("task_id")

    Task.configure(task_id, {"method": method})

    status_msg = Task.submit(task_id)

    return status_msg


@require_valid_task_id
@normalization_bp.route("/get-normalized-counts")
def get_normalized_counts():
    """Return the normalized counts file to the client"""

    #rel_user_dir = common.Task.get_dir(task_id)
    #abs_user_dir = os.path.abspath(rel_user_dir)
    #return send_from_directory(abs_user_dir, "counts_normalized.tsv")
    return "Not implemented", 404
