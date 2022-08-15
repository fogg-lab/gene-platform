import os
import subprocess
from ..common import common_blueprint as common
import time
from . import preprocessing_utils as utils
from flask import Blueprint, render_template, request, send_from_directory

# Set the current working directory and relative path to the user files
SCRIPT_PATH = os.path.realpath(__file__)
SCRIPT_DIR = "/".join(SCRIPT_PATH.split("/")[:-1])
USER_FILES_LOCATION = "../user_files"
PREP_GDC_SCRIPT = "Rscript ../../rscripts/prep_gdc.r"
PREP_GEO_SCRIPT = "Rscript ../../rscripts/prep_geo.r"

os.chdir(SCRIPT_DIR)

preprocessing_bp = Blueprint('preprocessing_bp', __name__,
                             template_folder='../templates',
                             static_folder='../static')


@preprocessing_bp.route("/preprocessing")
def preprocessing():
    '''Load preprocessing page'''

    common.ensure_session_dir()

    cur_uploads, all_uploads = common.list_user_files()

    return render_template(
        "preprocessing.html",
        cur_uploads=cur_uploads, all_uploads=all_uploads, title="Preprocessing"
    )


@preprocessing_bp.route("/submit-preprocessing", methods=["POST"])
def submit_preprocessing():
    data_source = request.form.get("source")
    dsets = request.form.get("dsets")
    valid_dsets, _ = get_valid_invalid_dsets(dsets, data_source)

    dsets = ""
    for dset in valid_dsets:
        dsets += f"{dset} "
    dsets = dsets.strip()

    user_dir = common.get_session_dir()

    src = "geo" if "geo" in data_source.lower() else "gdc"

    if src == "geo":
        expected_unmapped_counts_paths = []
        for dset in valid_dsets:
            expected_unmapped_counts_paths.append(
                os.path.join(user_dir, f"{dset}_counts_unmapped.tsv"))

    expected_paths = []
    for dset in valid_dsets:
        dset = dset.lower()
        expected_counts_path = os.path.join(user_dir,
                                            f"{dset}_counts_processed.tsv")
        expected_coldata_path = os.path.join(user_dir,
                                            f"{dset}_coldata_processed.tsv")
        expected_paths.append(expected_counts_path)
        expected_paths.append(expected_coldata_path)

    # Delete any previous preprocessing results
    for file in os.listdir(user_dir):
        if file.endswith("processed.tsv"):
            os.remove(os.path.join(user_dir, file))

    log = os.path.join(user_dir, "log")

    if "gdc" in data_source.lower():
        subprocess.Popen([f"{PREP_GDC_SCRIPT} {user_dir}/ {dsets} "
                          f"1> {log} 2>& 1"], shell=True)
    elif "geo" in data_source.lower():
        subprocess.Popen([f"{PREP_GEO_SCRIPT} {user_dir}/ {dsets} "
                          f"1> {log} 2>& 1"], shell=True)

    status_msg = "Preprocessing complete."

    if src == "geo":
        # Map probes to gene symbols once the r script is finished getting the unmapped counts
        unmapped_counts_exists = False
        while not unmapped_counts_exists:
            unmapped_counts_exists = True
            for path in expected_unmapped_counts_paths:
                if not os.path.isfile(path):
                    unmapped_counts_exists = False
            if not unmapped_counts_exists:
                time.sleep(0.25)
        utils.prep_geo_counts(user_dir)

    is_output = False
    while not is_output:
        is_output = True
        for expected_path in expected_paths:
            if not os.path.isfile(expected_path):
                is_output = False
        if not is_output:
            time.sleep(0.25)

    print("\n\n\n\nDone, returning status message.\n\n\n\n")

    return status_msg


@preprocessing_bp.route("/confirm-preprocessing-submission", methods=["POST"])
def confirm_preprocessing_submission():
    '''Validate submitted datasets and display datasets to load before submission'''

    data_source = request.form.get("source")
    dsets = request.form.get("dsets")

    # get the analysis formula to display for the user
    confirmation_message = get_preprocessing_confirmation_msg(dsets, data_source)

    return confirmation_message


@preprocessing_bp.route("/get-preprocessed-data")
def get_preprocessed_data():
    user_dir = common.get_session_dir()
    processed_data = utils.zip_preprocessed_data(user_dir)
    zip_fname = processed_data.split("/")[-1]
    return send_from_directory(user_dir, zip_fname)


def get_preprocessing_confirmation_msg(dsets, source):
    '''Return html confirmation message with datasets to load'''

    err_msg = f"<p><b>Error:</b> No provided datasets were recognized from source: {source}.</p>"

    if not dsets:
        return err_msg
    
    valid_dsets, invalid_dsets = get_valid_invalid_dsets(dsets, source)
    if not valid_dsets:
        return err_msg

    status_msg = ""

    if len(invalid_dsets) > 0:
        status_msg += "<p><b>Warning:</b> The following datasets were not recognized:"
        status_msg += "</b></p>"
        for dset in invalid_dsets:
            status_msg += f"<p>{dset}</p>"
        status_msg += "<br>"

    status_msg += f"<p><b>You are requesting the following datasets from {source}:"
    status_msg += "</b></p>"
    for dset in valid_dsets:
        status_msg += f"<p>{dset.upper()}</p>"
    status_msg += "<br>"
    status_msg += "<p>Proceed?</p>"

    return status_msg


def get_valid_invalid_dsets(dsets, source):
    if " " in dsets and not "," in dsets:
        dsets = dsets.split(" ")
    elif "," in dsets:
        dsets = dsets.split(",")
    else:
        dsets = [dsets.lower()]
    for i in range(len(dsets)):
        dsets[i] = dsets[i].replace(" ", "").lower()

    if "gdc" in source.lower():
        valid_dsets = utils.get_valid_gdc_projects(dsets)
    elif "geo" in source.lower():
        valid_dsets = utils.get_valid_geo_accessions(dsets)

    invalid_dsets = list(set(dsets).difference(set(valid_dsets)))

    return valid_dsets, invalid_dsets
