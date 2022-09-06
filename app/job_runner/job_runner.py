import os
from redis import Redis
from app.exceptions import JobNotFoundError
from app.job_runner.job_utils import (analysis, batch_correction, correlation,
                                      normalization, preprocessing,
                                      expression_data_validation)

JOB_TYPE_UTIL_MODULES = {
    "analysis": analysis,
    "batch_correction": batch_correction,
    "correlation": correlation,
    "normalization": normalization,
    "preprocessing": preprocessing,
    "expression_data_validation": expression_data_validation
}

JOB_TYPE_INPUT_FNAMES = {
    "analysis": analysis.INPUT_FNAMES,
    "batch_correction": batch_correction.INPUT_FNAMES,
    "correlation": correlation.INPUT_FNAMES,
    "normalization": normalization.INPUT_FNAMES,
    "preprocessing": preprocessing.INPUT_FNAMES
}

def run_job(job_dir, job_type):
    """Run a job, return status message"""
    status_msg = ""
    if job_type in JOB_TYPE_UTIL_MODULES:
        status_msg = JOB_TYPE_UTIL_MODULES[job_type].start_job(job_dir)
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

def add_input_file(file_contents, job_dir, standard_fname, user_fname, job_type):
    """
    Saves uploaded file for a job and perform input validation
    Args:
        file_contents: bytes object
        standard_filename: string
        user_filename: string
    Returns:
        status: dict
    """
    if job_type not in JOB_TYPE_UTIL_MODULES:
        raise Exception("Invalid job type")
    save_path = os.path.join(job_dir, "input", standard_fname)
    if os.path.exists(save_path):
        os.remove(save_path)
    lines = file_contents.split(b'\n')
    with open(save_path, "w", encoding="UTF-8") as user_file:
        for line in lines:
            user_file.write(f"{line.decode('UTF-8')}\n")
    # Perform input validation
    status = JOB_TYPE_UTIL_MODULES[job_type].update_job(save_path)
    if "ERROR" in status:
        os.remove(save_path)
    else:
        redis_db = Redis()
        job_id = get_job_id_from_dir(job_dir)
        redis_db.hset(f"{job_id}_input_files", standard_fname, user_fname)
    return status

def list_input_files(job_dir):
    """List base file names of all input files for a job"""
    redis_db = Redis()
    job_id = get_job_id_from_dir(job_dir)
    input_files = {}
    if redis_db.exists(job_id):
        input_files = redis_db.hgetall(f"{job_id}_input_files")
    else:
        raise JobNotFoundError("Job not found, so input files cannot be listed")

    return input_files

def get_job_id_from_dir(job_dir):
    return job_dir.rstrip("/").split("/")[-1]
