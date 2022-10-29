from pathlib import Path
from flask import (Blueprint, render_template, request, redirect, url_for,
                   Response, jsonify, send_from_directory)

from app.models.task import Task
from app.blueprints.common import require_valid_task_id
from app import helper

analysis_bp = Blueprint('analysis_bp', __name__)


@analysis_bp.route("/setup")
def setup():
    """first input page (file uploads)"""

    task_id = request.args.get("task_id")

    if not task_id:
        task_id = Task.create("analysis")

    uploads = Task.list_input_files(task_id)

    return render_template("uploads_form.html", uploaded_input_files=uploads,
                           title="DGE Analysis")


@require_valid_task_id
@analysis_bp.route("/upload", methods=["POST"])
def upload():
    """Receive uploaded file for analysis task"""

    user_fname = request.args.get("user_filename")
    standard_fname = request.headers.get('X_FILENAME')
    task_id = request.args.get("task_id")

    result = Task.add_input_file(task_id, request.data, standard_fname, user_fname)

    return jsonify(result)


@require_valid_task_id
@analysis_bp.route("/analysis-parameters", methods=["GET"])
def analysis_parameters():
    """Loads the analysis parameter form"""

    task_id = request.args.get("task_id")

    config_params = Task.get_config(task_id)

    return render_template("parameters_form.html", params=config_params,
                           title="Analysis parameters")


@require_valid_task_id
@analysis_bp.route("/confirm-analysis-submission", methods=["POST"])
def confirm_analysis_submission():
    """validate input and display formula before submission"""

    data_type = request.form.get("data_type")
    task_id = request.args.get("task_id")

    # Generate config file from form parameters
    params = get_analysis_request_parameters(request.form, data_type)
    status = Task.configure(task_id, params)

    # Validate config and input files
    if len(status.get("errors", [])) == 0:
        status = Task.validate_task(task_id)

    if len(status.get("errors", [])) > 0:
        for error in status["errors"]:
            confirmation_message = f"<p>Error: {error}</p>"
        return confirmation_message

    if len(status.get("warnings", [])):
        for warning in status["warnings"]:
            confirmation_message += f"<p>Warning: {warning}</p>"

    analysis_formula = status.get("status", "")
    confirmation_message += f"<p>{analysis_formula}</p>"

    return confirmation_message


@require_valid_task_id
@analysis_bp.route("/submit", methods=["POST"])
def submit():
    """Submit input for analysis"""

    task_id = request.args.get("task_id")

    status_msg = Task.submit(task_id)

    return status_msg


@require_valid_task_id
@analysis_bp.route("/display")
def display_output():
    """
    display output in a table
    """

    task_id = request.args.get("task_id")

    # check here if output.tsv exists and errors.txt doesn't
    unfiltered_output_path = Task.get_output_filepath(task_id, "output.tsv")

    unfiltered_output = helper.get_tsv_rows(unfiltered_output_path)
    if unfiltered_output is None:
        return "Analysis output file not found.", 204
    elif len(unfiltered_output) == 0:
        return "Analysis output file is empty.", 204
    else:
        unfiltered_output = unfiltered_output[:1]

    filter_output_path = Task.get_output_filepath(task_id, "filter_output.tsv")
    if filter_output_path == "":
        filtered_data = None    # no filter_output.tsv file
    else:
        filtered_output = helper.get_tsv_rows(filter_output_path)
        filtered_data = filtered_output[1:]

    return render_template("results.html", cols = unfiltered_output[0],
         data=unfiltered_output, filtered_data=filtered_data)


@require_valid_task_id
@analysis_bp.route("/plots")
def plots():
    """Load analysis plots page"""

    task_id = request.args.get("task_id")

    task_output_filepaths = Task.get_all_output_filepaths(task_id)
    plot_filenames = [Path(fp).name for fp in task_output_filepaths if Path(fp).suffix == ".png"]

    for filename in plot_filenames:
        if "mean" in filename and "unfiltered" in filename:
            plot_filenames["unfiltered_mean_variance"] = filename
        elif "mean" in filename:
            plot_filenames["filtered_mean_variance"] = filename
        elif "unfilt" in filename:
            plot_filenames["unfiltered_volcano"] = filename
        elif "filt" in filename:
            plot_filenames["filtered_volcano"] = filename

    return render_template("plots.html", plot_filenames=plot_filenames)


@require_valid_task_id
@analysis_bp.route('/getplot/<filename>')
def get_plot(filename):
    """Get an image for a plot"""

    task_id = request.args.get("task_id")

    plot_filepath = Task.get_output_filepath(task_id, filename)
    if plot_filepath == "":
        return f"Plot '{filename}' not found.", 204

    plot_dir = Path(plot_filepath).parent

    return send_from_directory(plot_dir, filename)


@analysis_bp.route("/reset", methods=["GET"])
def reset():
    """Redirects to the homepage without passing it the current task id"""
    return redirect(url_for("common_bp.index"))


@analysis_bp.route("/get-unfiltered-tsv")
def get_unfiltered_tsv():
    """Download unfiltered output."""

    task_id = request.args.get("task_id")

    unfiltered_output_path = Task.get_output_filepath(task_id, "output.tsv")

    if unfiltered_output_path == "":
        return "Analysis output file not found.", 204

    with open(unfiltered_output_path, encoding="UTF-8") as unfiltered_output:
        return Response(
            unfiltered_output,
            mimetype='text/csv',
            headers={'Content-disposition':
                    'attachment; filename=output.tsv'})


@analysis_bp.route("/get-filtered-tsv")
def get_filtered_tsv():
    """Download filtered output."""

    task_id = request.args.get("task_id")

    filtered_output_path = Task.get_output_filepath(task_id, "filter_output.tsv")

    if filtered_output_path == "":
        return "Filtered analysis output file not found.", 204

    with open(filtered_output_path, encoding="UTF-8") as filtered_output:
        return Response(
            filtered_output,
            mimetype="text/csv",
            headers={"Content-disposition":
                    "attachment; filename=filter_output.tsv"})


def get_analysis_request_parameters(form, data_type):
    """returns request parameters from the parameters form"""

    request_parameters = {}

    # if analysis type is microarray, consider use_qual_weights
    if data_type != "rnaseq":
        # convert 'None'|'on' to True|False
        if form.get("use_qual_weights") is None:
            request_parameters["use_qual_weights"] = False
        else:
            request_parameters["use_qual_weights"] = True

    request_parameters["min_prop"] = form.get("min_prop")
    request_parameters["min_expr"] = form.get("min_expr")
    request_parameters["adj_method"] = form.get("adj_method")
    request_parameters["contrast_level"] = form.get("contrast_level")
    request_parameters["reference_level"] = form.get("reference_level")
    request_parameters["padj_thresh"] = form.get("padj_thresh")

    return request_parameters
