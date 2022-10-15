import os
import shutil
from redis import Redis
import yaml
from app.exceptions import InvalidJobInputFile

class JobRunner:
    """Base class to prepare and run user jobs."""
    def __init__(self, job_id, job_dir):
        self._job_id = job_id
        self._job_dir = job_dir
        self._job_type = None
        self._input_filenames = None

    def configure(self, config):
        """
        Validate and save a config file for a job then return status message
        Args:
            config (dict): config parameters
        Returns:
            status (dict): status message with list of errors/warnings in "errors" and "warnings"
        """

        status = self.validate_config(config)
        if len(status["errors"]) == 0:
            # save config parameters to config.yml file
            config_path = os.path.join(self._job_dir, "input", "config.yml")
            with open(config_path, "w", encoding="utf-8") as cfg_file:
                yaml.dump(config, cfg_file)

        return status

    def get_log_update(self, last_log_update_line_number):
        """Returns the contents of the log file for a job"""
        log_file_path = os.path.join(self._job_dir, ".log")
        log_content = ""
        if os.path.isfile(log_file_path):
            with open(log_file_path, "r", encoding="utf-8") as log_file:
                # Limit to 100000 characters
                full_log_content = log_file.read()[-100000:]
                if len(full_log_content) > last_log_update_line_number:
                    log_content = full_log_content[last_log_update_line_number:]
        return log_content

    def add_input_file(self, file_contents, standard_fname, user_fname):
        """
        Saves uploaded file for a job and perform input validation
        Args:
            file_contents: bytes object
            standard_filename: string
            user_filename: string
        Returns:
            status: dict
        """
        if standard_fname not in self._input_filenames:
            raise InvalidJobInputFile(
                "Failed to add input file: "
                f"Filename '{standard_fname}' is invalid for job type: {self._job_type}. "
                f"Input file must be one of {self._input_filenames}")
        save_path = os.path.join(self._job_dir, "input", standard_fname)
        if os.path.exists(save_path):
            os.remove(save_path)
        lines = file_contents.split(b'\n')
        with open(save_path, "w", encoding="UTF-8") as user_file:
            for line in lines:
                user_file.write(f"{line.decode('UTF-8')}\n")
        # Perform input validation
        status = self.update_job(save_path)
        if "errors" in status:
            os.remove(save_path)
        else:
            redis_db = Redis()
            redis_db.hset(f"{self._job_id}_input_files", standard_fname, user_fname)
        return status

    def list_input_files(self):
        """List base file names of all input files for a job"""
        redis_db = Redis()
        input_files = {}
        redis_hash_name = f"{self._job_id}_input_files"
        if redis_db.exists(redis_hash_name):
            input_files = redis_db.hgetall(redis_hash_name)
        return input_files

    def get_job_id_from_dir(self):
        """Get job id from job directory, which is the basename of the directory"""
        return self._job_dir.rstrip("/").split("/")[-1]

    def remove_job(self):
        """Remove job directory and associated redis data"""
        redis_db = Redis()
        redis_hash_name = f"{self._job_id}_input_files"
        if redis_db.exists(redis_hash_name):
            all_keyval_pairs = redis_db.hgetall(redis_hash_name)
            all_keys = list(all_keyval_pairs.keys())
            redis_db.hdel(redis_hash_name, *all_keys)
        shutil.rmtree(self._job_dir)
