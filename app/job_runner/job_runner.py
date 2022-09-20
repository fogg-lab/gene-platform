"""The job runner module is used to prepare and run jobs."""

import os
import shutil
from redis import Redis
import yaml
from app.exceptions import InvalidJobType, InvalidInputFile
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


def configure_job(job_dir, job_type, config):
    """
    Validate and save a config file for a job then return status message
    Args:
        job_dir (str): path to job directory
        job_type (str): type of job (ie, "analysis")
        config (dict): config parameters
    """

    status = JOB_UTIL_MODULES[job_type].validate_config(config)
    if len(status["errors"]) == 0:
        # save config parameters to config.yml file
        config_path = os.path.join(job_dir, "input", "config.yml")
        with open(config_path, "w", encoding="utf-8") as cfg_file:
            yaml.dump(config, cfg_file)

    return status

def run_job(job_dir, job_type):
    """Run a job, return status message"""
    if job_type not in JOB_UTIL_MODULES:
        raise InvalidJobType(f"Failed to run job: Job type '{job_type}' is invalid."
                             f"Must be one of {JOB_UTIL_MODULES}")
    status_msg = JOB_UTIL_MODULES[job_type].start_job(job_dir)
    return status_msg

def get_job_log_update(job_dir, last_log_update_line_number):
    """Returns the contents of the log file for a job"""
    log_file_path = os.path.join(job_dir, ".log")
    log_content = ""
    if os.path.isfile(log_file_path):
        with open(log_file_path, "r", encoding="utf-8") as log_file:
            # Limit to 100000 characters
            full_log_content = log_file.read()[-100000:]
            if len(full_log_content) > last_log_update_line_number:
                log_content = full_log_content[last_log_update_line_number:]
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
    if job_type not in JOB_UTIL_MODULES:
        raise InvalidJobType(
            f"Failed to add input file: Job type '{job_type}' is not valid. "
            f"Job type must be one of {JOB_UTIL_MODULES}")
    if standard_fname not in JOB_UTIL_MODULES[job_type].INPUT_FNAMES:
        raise InvalidInputFile(
            "Failed to add input file: "
            f"Filename '{standard_fname}' is invalid for job type: {job_type}. "
            f"Input file must be one of {JOB_UTIL_MODULES[job_type].INPUT_FNAMES}")
    save_path = os.path.join(job_dir, "input", standard_fname)
    if os.path.exists(save_path):
        os.remove(save_path)
    lines = file_contents.split(b'\n')
    with open(save_path, "w", encoding="UTF-8") as user_file:
        for line in lines:
            user_file.write(f"{line.decode('UTF-8')}\n")
    # Perform input validation
    status = JOB_UTIL_MODULES[job_type].update_job(save_path)
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
    redis_hash_name = f"{job_id}_input_files"
    if redis_db.exists(redis_hash_name):
        input_files = redis_db.hgetall(redis_hash_name)
    return input_files

def validate_submission(job_dir, job_type):
    """Validate input files for submission and return status"""
    status = JOB_UTIL_MODULES[job_type].validate_submission(job_dir)
    return status

def get_job_id_from_dir(job_dir):
    """Get job id from job directory, which is the basename of the directory"""
    return job_dir.rstrip("/").split("/")[-1]

def remove_job(job_dir):
    """Remove job directory and associated redis data"""
    redis_db = Redis()
    job_id = get_job_id_from_dir(job_dir)
    redis_hash_name = f"{job_id}_input_files"
    if redis_db.exists(redis_hash_name):
        all_keyval_pairs = redis_db.hgetall(redis_hash_name)
        all_keys = list(all_keyval_pairs.keys())
        redis_db.hdel(redis_hash_name, *all_keys)
    shutil.rmtree(job_dir)
