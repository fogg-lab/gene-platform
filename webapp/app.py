import os
import subprocess
import shutil
import csv
from datetime import timedelta
import tempfile
from werkzeug.utils import secure_filename
from flask import Flask, render_template, request, redirect, url_for, \
    send_from_directory, session, Response
from flask_session.__init__ import Session

app = Flask(__name__)

app.config['SESSION_PERMANENT'] = True
app.config['SESSION_TYPE'] = 'filesystem'
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(minutes=30)
app.config['DROPZONE_TIMEOUT'] = 120000 # timeout for uploads in milliseconds

Session(app)

# Ensure that the current working directory is the webapp directory
# Get the path to the webapp dir from the path of this script
SCRIPT_PATH = os.path.realpath(__file__)
SCRIPT_DIR = '/'.join(SCRIPT_PATH.split('/')[:-1])
os.chdir(SCRIPT_DIR)

RNA_SEQ_SCRIPT = 'Rscript ../analysis/rnaseq.R'
MICROARRAY_SCRIPT = 'Rscript ../analysis/microarray.R'
USER_FILES_LOCATION = 'user_files'

# backlog (in order)
# - add validation of parameters
# - add ability to handle error output if analysis script fails

@app.route('/')
def index():
    files = os.listdir(USER_FILES_LOCATION)
    session_dir = get_session_dir()
    min_expr, min_prop, padj_thresh, adj_method, condition, contrast_level, reference_level = "","","","","","",""
    if session_dir:
        config_path = session["user_session_dir"] + "/config.yml"
        if os.path.exists(config_path):
            params = parse_config(config_path)
            min_expr, min_prop, padj_thresh, adj_method, condition, contrast_level, reference_level = params["min_expr"], params["min_prop"], params["padj_thresh"], params["adj_method"], params["condition"], params["contrast_level"], params["reference_level"]
            # return render_template('input_form_prefilled.html', files=files, min_expr=min_expr, min_prop=min_prop, padj_thresh=padj_thresh, adj_method=adj_method, condition=condition, contrast_level=contrast_level, reference_level=reference_level)

    return render_template('input_form.html', use_qual_weights=True)

@app.route('/uploadcounts', methods=['POST'])
def upload_counts():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        if file_ext != '.tsv':
            return 'Invalid file extension', 400
        save_temp_file(uploaded_file, 'counts.tsv')
    return redirect(url_for('index'))

@app.route('/uploadcoldata', methods=['POST'])
def upload_coldata():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        if file_ext != '.tsv':
            return 'Invalid file extension', 400
        save_temp_file(uploaded_file, 'coldata.tsv')
    return redirect(url_for('index'))

@app.route('/uploadfilter', methods=['POST'])
def upload_filter():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        if file_ext != '.txt':
            return 'Invalid file extension', 400
        save_temp_file(uploaded_file, 'filter.txt')
    return redirect(url_for('index'))

@app.route('/uploadconfig', methods=['POST'])
def upload_config():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        if file_ext != '.yml' and file_ext != '.txt':
            return 'Invalid file extension', 400
        save_temp_file(uploaded_file, 'config.yml')

    config_dict = parse_config('config.yml')
    if type(config_dict) == str:
        print(err_msg)
        # TODO: reload page with error message here
    else:
        print("config file is valid")

    return redirect(url_for('prefilled_index'))

@app.route('/user_files/<filename>')
def upload(filename):
    return send_from_directory('user_files', filename)

@app.route('/getconfig', methods=['POST'])
def get_config():
    parameters = {}
    if 'user_session_dir' in session:
        config_filepath = '%sconfig.yml' %(session['user_session_dir'])
    if 'user_session_dir' in session and os.path.isfile(config_filepath):
        config = open(('%sconfig.yml' %(session['user_session_dir'])), 'r')
        parameters = get_config_parameters(config)
        parameters["config_exists"] = True
    else:
        parameters["config_exists"] = False
    return parameters

# submit()
# submits input for analysis
# generates the config file if the parameters are filled out
@app.route('/submit', methods=['POST'])
def submit():
    # get whether analysis is microarray or RNA-Seq
    data_type = request.form.get('data_type')
<<<<<<< Updated upstream

=======

    # get config parameters
    min_expr = request.form.get('min_expr')
    min_prop = request.form.get('min_prop')
    padj_thresh = request.form.get('padj_thresh')
    adj_method = request.form.get('adj_method')
    condition = request.form.get('condition')
    contrast_level = request.form.get('contrast_level')
    reference_level = request.form.get('reference_level')

>>>>>>> Stashed changes
    # request.form.get('use_qual_weights') returns either "None" or "on"
    # needs to be boolean True or False
    use_qual_weights = request.form.get('use_qual_weights')

    parameters = {}
    # if analysis type is RNASeq, use padj_thresh
    if data_type == 'RNA-Seq':
        parameters['padj_thresh'] = request.form.get('padj_thresh')
    parameters['min_prop'] = request.form.get('min_prop')
    parameters['min_expr'] = request.form.get('min_expr')
    parameters['adj_method'] = request.form.get('adj_method')
    parameters['condition'] = request.form.get('condition')
    parameters['contrast_level'] = request.form.get('contrast_level')
    parameters['reference_level'] = request.form.get('reference_level')
    if use_qual_weights is None:
        parameters['use_qual_weights'] = False
    else:
        parameters['use_qual_weights'] = True

    parameters_filled = True
    for parameter_val in parameters:
        if not parameter_val:
            parameters_filled = False

        # if the parameters are set, generate config
    if parameters_filled:
        generate_config(parameters)

    if data_type == "microarray":
        subprocess.Popen(['%s %s' \
            %(MICROARRAY_SCRIPT, session["session_id"])], shell=True)
    elif data_type == "RNA-Seq":
        subprocess.Popen(['%s %s' \
            %(RNA_SEQ_SCRIPT, session["session_id"])], shell=True)

    # wait for the output.tsv file to appear in the session directory,
    # then redirect to the results page
    # TODO: error handling
    remove_old_output()
    analysis_done = False
    while not analysis_done:
        path_to_output =  session['user_session_dir'] + 'output.tsv'

        # '2>/dev/null' suppresses the expected error output 'file not found'
        output = subprocess.Popen(['ls %s %s' %(path_to_output, '2>/dev/null')\
            ], stdout=subprocess.PIPE, shell=True).communicate()[0]

        # results of the ls are returned in bytes, ends with newline character
        analysis_done = output == str.encode(path_to_output) + b'\n'

    return redirect(url_for('display_output'))

@app.route('/display')
def display_output():
    # check here if output.tsv exists and errors.txt doesn't
    path_to_output =  session['user_session_dir'] + 'output.tsv'
    output = open(path_to_output)
    reader = csv.reader(output, delimiter='\t')
    rows = [[elem for elem in row] for row in reader]
    output.close()

    # Delete the user session directory, including the input and output files
    cleanup_session()
    return render_template('results.html', col_titles=rows[:1], info=rows[1:])

@app.route('/getunfilteredtsv')
def get_unfiltered_tsv():
    unfiltered_output = open("%soutput.tsv" %(session['user_session_dir']))
    return Response(
        unfiltered_output,
        mimetype="text/csv",
        headers={"Content-disposition":
                "attachment; filename=output.tsv"})

@app.route('/getfilteredtsv')
def get_filtered_tsv():
    filtered_output = open("%sfilter_output.tsv" %(session['user_session_dir']))
    return Response(
        filtered_output,
        mimetype="text/csv",
        headers={"Content-disposition":
                "attachment; filename=filter_output.tsv"})

# Takes a user's file and copies it into a temp directory on the server
# directory path is stored in the user session variable "user_session_dir"
def save_temp_file(file, filename):
    if ('user_session_dir' not in session or
            not os.path.exists(session['user_session_dir'])):
        temp_dir = tempfile.mkdtemp(dir=USER_FILES_LOCATION)
        os.chmod(temp_dir, 0o777) # give everyone rwx permission for the dir
        session['user_session_dir'] = temp_dir + '/'
        session['session_id'] = temp_dir.split('/')[-1:]
        schedule_session_deletion(session['user_session_dir'])

    user_file_path = session['user_session_dir'] + filename

    if os.path.exists(user_file_path):
        os.remove(user_file_path)
    user_file = open(user_file_path, 'w')

    # copy the file line by line, adding a trailing newline if none exists
    lines = file.readlines()
    for line in lines:
        user_file.write(line.decode('utf-8'))

    user_file.close()

# cleanup_session() cleans up the user's session
# Removes the temp directory for the session, including input/output files
def cleanup_session():
    if 'user_session_dir' in session:
        pass
        # need to find a different way, this deletes output files prematurely
        #shutil.rmtree(session['user_session_dir'])

def cleanup_old_sessions():
    # Clean up other sessions older than one day (86400 seconds)
    for old_dir in os.listdir(USER_FILES_LOCATION):
        old_dir = '%s/%s' %(USER_FILES_LOCATION, old_dir)

        if old_dir != '%s/.gitkeep' %(USER_FILES_LOCATION):
            get_age = '$(($(date +%s) - $(date +%s -r ' + old_dir + ')))'

            # Run bash command and return the stdout output
            age_seconds = int(subprocess.Popen(['echo %s' %(get_age)], \
                stdout=subprocess.PIPE, shell=True).communicate()[0])

            if age_seconds > 86400:
                shutil.rmtree(old_dir)

def get_session_dir():
    if "user_session_dir" in session:
        return session["user_session_dir"]
    else:
        return False

def generate_config(parameters):

    config_file_path = session["user_session_dir"] + "/config.yml"
    config_file = open(config_file_path, 'w')

    for param in parameters.keys():
        if param in ['adj_method', 'condition', 'contrast_level',\
             'reference_level']:
            config_file.write('%s: \"%s\"\n' %(param, parameters[param]))
        else:
            config_file.write('%s: %s\n' %(param, parameters[param]))

    config_file.close()

def get_config_parameters(config_file):
    parameters = {}
    config_lines = config_file.readlines()
    for line in config_lines:
        # remove quotes and newline characters
        line = line.translate(str.maketrans('', '', '\n\"\''))
        if line:
            parameter_name, parameter_value = tuple(line.split(': '))
            parameters[parameter_name] = parameter_value
    return parameters

# removes old output files from the current session
def remove_old_output():
    filtered_output_path = "%sfilter_output.tsv" %(session['user_session_dir'])
    unfiltered_output_path = "%soutput.tsv" %(session['user_session_dir'])
    filtered_output_exists = os.path.exists(filtered_output_path)
    unfiltered_output_exists = os.path.exists(unfiltered_output_path)
    if filtered_output_exists:
        os.remove(filtered_output_path)
    if unfiltered_output_exists:
        os.remove(unfiltered_output_path)

def schedule_session_deletion(session_path):
    pass

def parse_config(filename):
    params = {}
    config_file_path = session["user_session_dir"] + "/config.yml"
    file_text = (open(config_file_path, 'r')).readlines()
    # params = yaml.load(file_text)
    for line in file_text:
        # try:
        keyvalue = line.split(': ')
        key = keyvalue[0]
        value = keyvalue[1][:-1] #remove trailing newline
        if is_numeric(value):
            value = float(value)
        elif is_bool(value):
             value = bool(value)
        else:
            value = value[1:-1]
        params[key] = value
        # except:
        #     return "config file is not formatted correctly"

    err_msg = validate_config_params(params)
    if err_msg != "":
        return err_msg
    return params

def is_numeric(string):
    num_symbols = "-.0123456789"
    return all([char in num_symbols for char in string])

def is_bool(string):
    return string in {"True", "False"}

# returns an error message if config parameters are invalid
# otherwise, returns an empty string
def validate_config_params(params):
    if type(params["min_expr"]) not in [int, float]:
        return "min_expr must be a number"
    if params["min_expr"] < 0:
        return "min_expr must be a non-negative"
    if type(params["min_prop"]) not in [int, float]:
        return "min_prop must be a number"
    if params["min_prop"] < 0:
        return "min_prop must be a non-negative"
    if type(params["padj_thresh"]) not in [int, float]:
        return "padj_thresh must be a number"
    if type(params["adj_method"]) != str:
        return "adj_method must be a string"
    if type(params["condition"]) != str:
        return "condition must be a string"
    if type(params["contrast_level"]) != str:
        return "contrast_level must be a string"
    if type(params["reference_level"]) != str:
        return "reference_level must be a string"
    if type(params["use_qual_weights"]) != bool:
        return "reference_level must be a bool"
    return ""

cleanup_old_sessions()
