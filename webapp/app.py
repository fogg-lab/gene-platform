import os
import csv
from flask import Flask, render_template, request, redirect, url_for, \
    send_from_directory
from werkzeug.utils import secure_filename

app = Flask(__name__)

# backlog (in order)
# - limit upload fields to one file (Noah)
# - add uploading of config file (Noah)
# - add display page html (Corbin)
# - add route function which parses basic .tsv output file and renders display page (Carter)
# - add user session implementation for uploading files and displaying output
# - add parameter text input fields
# - add ability to generate config file in correct session directory
# - add submit button which generates config and calls analysis script
# - add validation of parameters and uploaded files
# - add ability to handle error.tsv output if analysis script fails
# - add ability to choose microarray or RNA-seq
# - add default values to text input field

def validate_counts(stream):
    # Stub function for demo
    return '.tsv'

def validate_coldata(stream):
    # Stub function for demo
    return '.tsv'

def validate_filter(stream):
    # Stub function for demo
    return '.txt'

@app.route('/')
def index():
    files = os.listdir('uploads')
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
        uploaded_file.save('uploads/counts.tsv')
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
        uploaded_file.save('uploads/coldata.tsv')
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
        uploaded_file.save('uploads/filter.txt')
    return redirect(url_for('index'))

@app.route('/uploads/<filename>')
def upload(filename):
    return send_from_directory('uploads', filename)

@app.route('/display')
def display_output():
    # check here if output.tsv exists and errors.txt doesn't
    output = open('output.tsv')
    reader = csv.reader(output, delimiter='\t')
    rows = [[elem for elem in row] for row in reader]
    output.close()
    return render_template('results.html', rows=rows)
