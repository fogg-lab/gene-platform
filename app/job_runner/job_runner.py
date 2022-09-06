import os
from app.job_runner.job_utils import (analysis, batch_correction, correlation,
                                      normalization, preprocessing,
                                      expression_data_validation)

JOB_UTIL_MODULES = {
    "analysis": analysis,
    "batch_correction": batch_correction,
    "correlation": correlation,
    "normalization": normalization,
    "preprocessing": preprocessing,
    "expression_data_validation": expression_data_validation
}

def run_job(job_dir, job_type):
    """Run a job, return status message"""
    status_msg = ""
    if job_type in JOB_UTIL_MODULES:
        status_msg = JOB_UTIL_MODULES[job_type].start_job(job_dir)
    else:
        raise Exception("Invalid job type")
    return status_msg

def get_job_log_update(job_dir):
    """
    returns the contents of the log file in user session directory
    the log file contains terminal output from the analysis script
    """
    log_file_path = os.path.join(job_dir, ".log")
    log_content = ""
    if os.path.isfile(log_file_path):
        with open(log_file_path, "r", encoding="utf-8") as log_file:
            log_content = log_file.read()[-1000:]
        with open(log_file_path, "r+", encoding="utf-8") as log_file:
            log_file.truncate(0)
    return log_content

def add_input_file(file_contents, save_path, job_type):
    """
    Saves uploaded file for a job and perform input validation
    Args:
        file_contents: bytes object
        standard_filename: string
        user_filename: string
    """
    if job_type not in JOB_UTIL_MODULES:
        raise Exception("Invalid job type")
    if os.path.exists(save_path):
        os.remove(save_path)
    lines = file_contents.split(b'\n')
    with open(save_path, "w", encoding="UTF-8") as user_file:
        for line in lines:
            user_file.write(f"{line.decode('UTF-8')}\n")
    # Perform input validation
    status_msg = JOB_UTIL_MODULES[job_type].update_job(save_path)
    return status_msg

def list_input_files(job_dir):
    """List base file names of all input files for a job"""
    input_dir = os.path.join(job_dir, "input")
    input_files = [x.split(".")[0] for x in os.listdir(input_dir)]
    return input_files
