import os
import subprocess
import shutil
import csv
from datetime import timedelta
import tempfile
from werkzeug.utils import secure_filename
from flask import Flask, render_template, request, redirect, url_for, \
    send_from_directory, session
from flask_session.__init__ import Session

app = Flask(__name__)

app.config['SESSION_PERMANENT'] = True
app.config['SESSION_TYPE'] = 'filesystem'
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(minutes=30)
app.config['DROPZONE_TIMEOUT'] = 120000 # timeout for uploads in milliseconds
Session(app)

RNA_SEQ_SCRIPT_LOCATION = "dge_analysis/rna_seq_dgea.R"
MICROARRAY_SCRIPT_LOCATION = "dge_analysis/microarray_dgea.R"

# backlog (in order)
# - limit upload fields to one file (Noah) - done (a new upload overwrites old)
# - add uploading of config file (Noah) - done
# - add display page html (Corbin)
# - add route function which parses basic .tsv output file and renders
#    display page (Carter)
# - add user session implementation for uploading files and displaying output
# - add parameter text input fields
# - add ability to generate config file in correct session directory
# - add submit button which generates config and calls analysis script
# - add validation of parameters and uploaded files
# - add ability to handle error output if analysis script fails
# - add ability to choose microarray or RNA-seq

@app.route('/')
def index():
    files = os.listdir('user_files')
    return render_template('input_form.html', files=files)

@app.route('/uploadcounts', methods=['POST'])
def upload_counts():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        if file_ext != '.tsv':
            return "Invalid file extension", 400
        save_temp_file(uploaded_file, 'counts.tsv')
    return redirect(url_for('index'))

@app.route('/uploadcoldata', methods=['POST'])
def upload_coldata():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        if file_ext != '.tsv':
            return "Invalid file extension", 400
        save_temp_file(uploaded_file, 'coldata.tsv')
    return redirect(url_for('index'))

@app.route('/uploadfilter', methods=['POST'])
def upload_filter():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        if file_ext != '.txt':
            return "Invalid file extension", 400
        save_temp_file(uploaded_file, 'filter.txt')
    return redirect(url_for('index'))

@app.route('/uploadconfig', methods=['POST'])
def upload_config():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        if file_ext != '.txt' or file_ext != '.yml':
            return "Invalid file extension", 400
        save_temp_file(uploaded_file, 'config.yml')
    return redirect(url_for('index'))

@app.route('/user_files/<filename>')
def upload(filename):
    return send_from_directory('user_files', filename)

# submit()
# submits input for analysis
# generates the config file if the parameters are filled out
@app.route('/submit', methods=['POST'])
def submit():
    # get whether analysis is microarray or RNA-Seq
    data_type = request.args.get('data_type')
    
    # get config parameters
    min_expr = request.args.get('min_expr')
    min_prop = request.args.get('min_prop')
    padj_thresh = request.args.get('padj_thresh')
    adj_method = request.args.get('adj_method')
    condition = request.args.get('condition')
    contrast_level = request.args.get('contrast_level')
    reference_level = request.args.get('reference_level')
    use_qual_weights = request.args.get('use_qual_weights')

    # if the parameters are set, generate config
    if min_expr is not None and min_prop is not None \
        and padj_thresh is not None and adj_method is not None \
        and adj_method is not None and condition is not None \
        and contrast_level is not None and reference_level is not None \
        and use_qual_weights is not None:
        generate_config(min_expr, min_prop, padj_thresh, adj_method, condition, \
            contrast_level, reference_level, use_qual_weights)

    # input files should be in the user directory by now
    # TODO: validate input files before calling analysis

    if data_type == "microarray":
        subprocess.Popen(['MICROARRAY_SCRIPT_LOCATION %s' \
            %("webapp/" + session["user_session_dir"])], shell=True)
    elif data_type == "RNA_Seq":
        subprocess.Popen(['RNA_SEQ_SCRIPT_LOCATION %s' \
            %("webapp/" + session["user_session_dir"])], shell=True)

    return

@app.route('/display')
def display_output():
    # check here if output.tsv exists and errors.txt doesn't
    output = open('output.tsv')
    reader = csv.reader(output, delimiter='\t')
    rows = [[elem for elem in row] for row in reader]
    output.close()
    cleanup_session()
    return render_template('results.html', rows=rows)

# Takes a user's file and copies it into a temp directory on the server
# directory path is stored in the user session variable "user_session_dir"
def save_temp_file(file, filename):
    if ("user_session_dir" not in session or
            not os.path.exists(session["user_session_dir"])):
        temp_dir = tempfile.mkdtemp(dir="user_files")
        session["user_session_dir"] = temp_dir + '/'

    user_file_path = session["user_session_dir"] + filename

    if os.path.exists(user_file_path):
        os.remove(user_file_path)
    user_file = open(user_file_path, 'w')

    for line in file.readlines():
        user_file.write(line.decode("utf-8"))

    user_file.close()

# cleanup_session() cleans up the user's session on exit
# Removes the temp directory for the session, including input/output files
def cleanup_session():
    if "user_session_dir" in session:
        shutil.rmtree(session["user_session_dir"])

def cleanup_old_sessions():
    # Clean up other sessions older than one hour (3600 seconds)
    for old_dir in os.listdir('user_files'):
        old_dir = 'user_files/' + old_dir

        if old_dir != "user_files/.gitkeep":
            get_age = "$(($(date +%s) - $(date +%s -r " + old_dir + ")))"

            # Run bash command and return the stdout output
            age_seconds = int(subprocess.Popen(['echo %s' %(get_age)], \
                stdout=subprocess.PIPE, shell=True).communicate()[0])

            if age_seconds > 3600:
                shutil.rmtree(old_dir)

def generate_config(min_expr, min_prop, padj_thresh, adj_method, condition, \
        contrast_level, reference_level, use_qual_weights):
    config_file_path = session["user_session_dir"] + "/config.yml"
    config_file = open(config_file_path, 'w')

    config_file.write("min_expr:", str(min_expr))
    config_file.write("min_prop:", str(min_prop))
    config_file.write("padj_thresh:", str(padj_thresh))
    config_file.write("adj_method:", str(adj_method))
    config_file.write("condition:", str(condition))
    config_file.write("contrast_level:", str(contrast_level))
    config_file.write("reference_level:", str(reference_level))
    config_file.write("use_qual_weights:", str(use_qual_weights))


cleanup_old_sessions()
