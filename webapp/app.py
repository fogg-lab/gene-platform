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
    return render_template('uploadfiles.html', files=files)

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
        if file_ext != '.txt':
            return "Invalid file extension", 400
        save_temp_file(uploaded_file, 'config.txt')
    return redirect(url_for('index'))

@app.route('/user_files/<filename>')
def upload(filename):
    return send_from_directory('user_files', filename)

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
        user_file.write(line.decode("utf-8") + '\n')

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

cleanup_old_sessions()
