import os
from flask import Flask, render_template, request, redirect, url_for, abort, \
    send_from_directory
from werkzeug.utils import secure_filename

app = Flask(__name__)

# TODO: Only allow 1 file per upload box
# TODO: Store uploaded files separately for each user session
# TODO: Add ability to upload config file instead of entering parameters
# TODO: Add form to enter parameters, on submit it generates the config file
# TODO: Validate config file & other uploaded files

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
