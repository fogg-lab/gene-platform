import os
import csv
from flask import (Blueprint, render_template, request, redirect, url_for,
                   Response, jsonify, send_from_directory)
from app.models.job import Job
from app.job_runner import job_runner
from app import helper

analysis_bp = Blueprint('analysis_bp', __name__)

RNASEQ_SCRIPT = "dge_rnaseq.r"
MICROARRAY_SCRIPT = "dge_microarray.r"
JOB_TYPE = "analysis"

@analysis_bp.route("/setup")
def setup():
    """first input page (file uploads)"""

    job_id = request.args.get("job_id")

    if not job_id:
        job_id = Job.create(JOB_TYPE)
    job_dir = Job.get_dir(job_id)
    uploads = job_runner.list_input_files(job_dir)

    return render_template("uploads_form.html", cur_uploads=uploads,
                           title="DGE Analysis")


@analysis_bp.route("/upload", methods=["POST"])
def upload():
    """Receive uploaded file for analysis job"""

    result = {}

    user_fname = request.args.get("user_filename")
    standard_fname = request.headers.get('X_FILENAME')
    job_id = request.args.get("job_id")

    if not job_id:
        return redirect(url_for("analysis_bp.setup"))
    job_dir = Job.get_dir(job_id)

    result = job_runner.add_input_file(request.data, job_dir, standard_fname,
                                       user_fname, JOB_TYPE)

    return jsonify(result)


@analysis_bp.route("/analysis-parameters", methods=["GET"])
def analysis_parameters():
    """Loads the analysis parameter form"""

    job_id = request.args.get("job_id")

    if not job_id:
        return redirect(url_for("analysis_bp.setup"))
    job_dir = Job.get_dir(job_id)

    params = helper.yaml_to_dict(os.path.join(job_dir, "config.yml"))
    return render_template("parameters_form.html", params=params,
                           title="Analysis parameters")


@analysis_bp.route("/submit", methods=["POST"])
def submit():
    """Submit input for analysis"""

    job_id = request.args.get("job_id")

    if not job_id:
        return redirect(url_for("analysis_bp.setup"))
    result = Job.submit_job(job_id)

    return result


@analysis_bp.route("/confirm-analysis-submission", methods=["POST"])
def confirm_analysis_submission():
    """validate input and display formula before submission"""

    validation_status_messages = []

    data_type = request.form.get("data_type")
    job_id = request.args.get("job_id")
    if not job_id:
        return redirect(url_for("analysis_bp.setup"))
    job_dir = Job.get_dir(job_id)
    
    # Generate config file from form parameters
    params = get_analysis_request_parameters(request.form, data_type)
    params_yaml = helper.dict_to_yaml(params)
    config_status = job_runner.add_input_file(params_yaml, job_dir, "config.yml",
                                   "config.yml", JOB_TYPE)
    validation_status_messages.append(config_status)

    # Validate files for submission
    validation_status_messages += job_runner.validate_submission(job_dir, JOB_TYPE)

    # get the analysis formula to display for the user
    confirmation_message = ""

    for msg in validation_status_messages:
        confirmation_message += f"<p>{msg}</p>"

    return confirmation_message


@analysis_bp.route("/display")
def display_output():
    """
    display output in a table
    """

    job_id = request.args.get("job_id")
    if not job_id:
        return redirect(url_for("analysis_bp.setup"))
    job_dir = Job.get_dir(job_id)

    # check here if output.tsv exists and errors.txt doesn't
    unfiltered_output_path = os.path.join(job_dir, "output.tsv")
    unfiltered_output_file = open(unfiltered_output_path, encoding="UTF-8")
    output_reader = csv.reader(unfiltered_output_file, delimiter="\t")
    unfiltered_output = list(output_reader)
    unfiltered_output_file.close()

    filter_output_path = os.path.join(job_dir, "filter_output.tsv")
    filter_output_exists = os.path.exists(filter_output_path)
    filtered_data = None
    if filter_output_exists:
        filter_output_file = open(filter_output_path, encoding="UTF-8")
        output_reader = csv.reader(filter_output_file, delimiter="\t")
        filtered_output = list(output_reader)
        filter_output_file.close()
        filtered_data = filtered_output[1:]

    return render_template("results.html", cols = unfiltered_output[:1][0],
         data=unfiltered_output[1:], filtered_data=filtered_data)


@analysis_bp.route("/plots")
def plots():
    """Load analysis plots page"""

    job_id = request.args.get("job_id")
    if not job_id:
        return redirect(url_for("analysis_bp.setup"))
    job_dir = Job.get_dir(job_id)
    output_dir = os.path.join(job_dir, "output")

    plot_filenames = {}

    for filename in os.listdir(output_dir):
        if "mean" in filename and ".png" and "unfiltered" in filename:
            plot_filenames["unfiltered_mean_variance"] = filename
        elif "mean" in filename and ".png" in filename:
            plot_filenames["filtered_mean_variance"] = filename
        elif "unfilt" in filename and ".png" in filename:
            plot_filenames["unfiltered_volcano"] = filename
        elif ".png" in filename:
            plot_filenames["filtered_volcano"] = filename

    return render_template("plots.html", plot_filenames=plot_filenames)


@analysis_bp.route('/getplot/<filename>')
def get_plot(filename):
    """Get an image for a plot"""

    job_id = request.args.get("job_id")
    if not job_id:
        return redirect(url_for("analysis_bp.setup"))
    job_dir = Job.get_dir(job_id)
    output_dir = os.path.join(job_dir, "output")

    return send_from_directory(output_dir, filename)


@analysis_bp.route("/reset", methods=["GET"])
def reset():
    """Deletes current job and redirects to homepage"""

    job_id = request.args.get("job_id")
    if job_id:
        job_dir = Job.get_dir(job_id)
        Job.delete(job_id)
        job_runner.remove_job(job_dir)

    return redirect(url_for("common_bp.index"))


@analysis_bp.route("/get-unfiltered-tsv")
def get_unfiltered_tsv():
    """download unfiltered output"""

    job_id = request.args.get("job_id")
    if not job_id:
        return redirect(url_for("analysis_bp.setup"))
    job_dir = Job.get_dir(job_id)
    output_dir = os.path.join(job_dir, "output")

    unfiltered_output_path = os.path.join(output_dir, "output.tsv")
    with open(unfiltered_output_path, encoding="UTF-8") as unfiltered_output:
        return Response(
            unfiltered_output,
            mimetype='text/csv',
            headers={'Content-disposition':
                    'attachment; filename=output.tsv'})


@analysis_bp.route("/get-filtered-tsv")
def get_filtered_tsv():
    """download filtered output"""

    job_id = request.args.get("job_id")
    if not job_id:
        return redirect(url_for("analysis_bp.setup"))
    job_dir = Job.get_dir(job_id)
    output_dir = os.path.join(job_dir, "output")

    filtered_output_path = os.path.join(output_dir, "filter_output.tsv")
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
        # form.get("use_qual_weights") will initially be either 'None' or 'on'
        # it needs to be a boolean True or False
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
