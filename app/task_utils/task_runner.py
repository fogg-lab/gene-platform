import os
import shutil
from abc import ABC, abstractmethod
from redis import Redis
import yaml

from app.exceptions import InvalidTaskInputFile


class TaskRunner(ABC):
    """Abstract base class to prepare and run user tasks."""
    def __init__(self, task_id, task_dir):
        self._task_id = task_id
        self._task_dir = task_dir
        self._task_type = None
        self._input_filenames = None

    @abstractmethod
    def update_task(self):
        pass

    @abstractmethod
    def start_task(self):
        pass

    @abstractmethod
    def validate_config(self, config):
        pass

    def save_config(self, config):
        # save config parameters to config.yml file
        config_path = os.path.join(self._task_dir, "input", "config.yml")
        with open(config_path, "w", encoding="utf-8") as cfg_file:
            yaml.dump(config, cfg_file)

    def get_log_update(self, last_log_update_line_number):
        """Returns the contents of the log file for a task"""
        log_file_path = os.path.join(self._task_dir, ".log")
        log_content = ""
        if os.path.isfile(log_file_path):
            with open(log_file_path, "r", encoding="utf-8") as log_file:
                # Limit to 100000 characters
                full_log_content = log_file.read()[-100000:]
                if len(full_log_content) > last_log_update_line_number:
                    log_content = full_log_content[last_log_update_line_number:]
        return log_content

    def add_input_file(self, file_contents, standard_fname):
        """Saves uploaded file for a task and returns the save path."""

        if standard_fname not in self._input_filenames:
            raise InvalidTaskInputFile(
                "Failed to add input file: "
                f"Filename '{standard_fname}' is invalid for task type: {self._task_type}. "
                f"Input file must be one of {self._input_filenames}")

        save_path = os.path.join(self._task_dir, "input", standard_fname)

        if os.path.exists(save_path):
            os.remove(save_path)

        lines = file_contents.split(b'\n')

        with open(save_path, "w", encoding="UTF-8") as user_file:
            for line in lines:
                user_file.write(f"{line.decode('UTF-8')}\n")

        return save_path

    def list_input_files(self):
        """
        List base filenames of all input files for a task.
        Returns:
            list[str]: List of base filenames.
        """
        redis_db = Redis()
        input_files = {}
        redis_hash_name = f"{self._task_id}_input_files"
        if redis_db.exists(redis_hash_name):
            input_files = redis_db.hgetall(redis_hash_name)
        return input_files

    def remove_task(self):
        """Remove task directory and associated redis data"""
        redis_db = Redis()
        redis_hash_name = f"{self._task_id}_input_files"
        if redis_db.exists(redis_hash_name):
            all_keyval_pairs = redis_db.hgetall(redis_hash_name)
            all_keys = list(all_keyval_pairs.keys())
            redis_db.hdel(redis_hash_name, *all_keys)
        shutil.rmtree(self._task_dir)
