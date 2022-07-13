import os
from ..common import common_blueprint as common
from flask import Blueprint, render_template

# Set the current working directory and relative path to the user files
SCRIPT_PATH = os.path.realpath(__file__)
SCRIPT_DIR = "/".join(SCRIPT_PATH.split("/")[:-1])
USER_FILES_LOCATION = "../user_files"

os.chdir(SCRIPT_DIR)

preprocessing_bp = Blueprint('preprocessing_bp', __name__, template_folder='../templates', \
    static_folder='../static')

@preprocessing_bp.route("/preprocessing")
def preprocessing():
    '''todo'''

    common.ensure_session_dir()

    cur_uploads, all_uploads = common.list_user_files()

    return render_template("preprocessing.html", \
        cur_uploads=cur_uploads, all_uploads=all_uploads, title="Batch Correction")
