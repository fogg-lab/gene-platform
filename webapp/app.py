import os
import csv
from flask import Flask, render_template, request, redirect, url_for, \
    send_from_directory, session
from flask_session.__init__ import Session
from werkzeug.utils import secure_filename
from datetime import timedelta
import tempfile

app = Flask(__name__)

app.config['SESSION_PERMANENT'] = True
app.config['SESSION_TYPE'] = 'filesystem'
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(minutes=30)
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
# - add ability to handle error.tsv output if analysis script fails
# - add ability to choose microarray or RNA-seq
# - add default values to text input field

def validate_counts(stream):
    # TODO: implement validate_counts
    return '.tsv'

def validate_coldata(stream):
    # TODO: implement validate_coldata
    return '.tsv'

def validate_filter(stream):
    # TODO: implement validate_filter
    return '.txt'

def validate_config(stream):
    # TODO: implement validate_filter
    return '.txt'

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
        if file_ext != '.tsv' or \
                file_ext != validate_counts(uploaded_file.stream):
            return "Invalid file extension", 400
        save_temp_file(uploaded_file, 'counts')
    return redirect(url_for('index'))

@app.route('/uploadcoldata', methods=['POST'])
def upload_coldata():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        if file_ext != '.tsv' or \
                file_ext != validate_coldata(uploaded_file.stream):
            return "Invalid file extension", 400
        save_temp_file(uploaded_file, 'coldata')
    return redirect(url_for('index'))

@app.route('/uploadfilter', methods=['POST'])
def upload_filter():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        if file_ext != '.txt' or \
                file_ext != validate_filter(uploaded_file.stream):
            return "Invalid file extension", 400
        save_temp_file(uploaded_file, 'filter')
    return redirect(url_for('index'))

@app.route('/uploadconfig', methods=['POST'])
def upload_config():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        if file_ext != '.txt' or \
                file_ext != validate_config(uploaded_file.stream):
            return "Invalid file extension", 400
        save_temp_file(uploaded_file, 'config')
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
    return render_template('results.html', rows=rows)

# Takes a user's file and copies it into a temp directory on the server
# directory path is stored in the user session variable "user_files_dirname"
def save_temp_file(file, filename):
    if not "user_files_dirname" in session:
        temp_dirname = tempfile.TemporaryDirectory(dir='user_files/')
        session["user_files_dirname"] = temp_dirname.name + '/'
    user_file_path = session["user_files_dirname"] + filename
    
    if os.path.exists(user_file_path):
        os.remove(user_file_path)
    user_file = open(user_file_path, 'wb')

    for line in file:
        user_file.write(line)
