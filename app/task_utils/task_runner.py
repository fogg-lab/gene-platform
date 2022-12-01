import os
import shutil
from pathlib import Path
import subprocess
from itertools import groupby
from abc import ABC, abstractmethod
from typing import List, Tuple, Dict
from glob import glob
from redis import Redis
import yaml

from app.exceptions import InvalidTaskInputFile


class TaskRunner(ABC):
    """Abstract base class to prepare and run user tasks."""
    def __init__(self, task_id, task_dir):
        self._task_id = task_id
        self._task_dir = task_dir
        self.task_type = None
        self._input_filenames = None
        self._full_log_content = ""

    @abstractmethod
    def execute_task(self) -> dict:
        pass

    @abstractmethod
    def validate_config(self, config) -> dict:
        """Validates config file (implemented in child classes)."""
        pass

    @abstractmethod
    def validate_task(self) -> dict:
        """Validates task input - config and uploaded files (implemented in child classes)."""
        pass

    def get_output_filepath(self, standard_filename) -> str:
        """
        Gets the filepath of an output file for a task.
        Args:
            standard_filename (str): The standard filename of the output file.
        Returns:
            str: The filepath of the output file, or empty string if file not found.
        """
        out_path = os.path.join(self._task_dir, "output", standard_filename)

        if not os.path.isfile(out_path):
            out_path = ""

        return out_path

    def get_all_output_filepaths(self) -> List[str]:
        """
        Gets the filepaths of all output files for a task.
        Returns:
            List[str]: List of filepaths.
        """
        out_paths = [out_path for out_path in glob(os.path.join(self._task_dir, "output", "*"))]
        out_paths = [fp for fp in out_paths if os.path.isfile(fp) and not fp.endswith(".zip")]
        return out_paths

    def get_config(self):
        """Returns the config parameters for the task as a dict."""
        config_path = os.path.join(self._task_dir, "input", "config.yml")
        if not os.path.isfile(config_path):
            return {}
        with open(config_path, "r", encoding="utf-8") as cfg_file:
            return yaml.safe_load(cfg_file)

    def save_config(self, config):
        """Saves config parameters to a config.yml file."""
        config_path = os.path.join(self._task_dir, "input", "config.yml")
        with open(config_path, "w", encoding="utf-8") as cfg_file:
            yaml.dump(config, cfg_file)

    @staticmethod
    def _filter_log(log_path, task_type):
        """Filters log content to remove extraneous output."""

        with open(log_path, "r") as log_file:
            full_log_content = log_file.read()

        # Check first 20 lines for extraneous warnings and remove them
        extraneous_output_substrings = [
            "During startup - Warning messages",
            "Setting LC_",
            "System has not been booted with systemd as init system",
            "Failed to create bus connection",
            "[1] TRUE",
            "null device",
        ]
        extraneous_line_indices = []

        log_lines = full_log_content.splitlines()
        for line_idx, line in enumerate(log_lines[:20]):
            for substring in extraneous_output_substrings:
                if substring in line:
                    extraneous_line_indices.append(line_idx)

        for line_idx in sorted(extraneous_line_indices, reverse=True):
            del log_lines[line_idx]

        # For consecutive download progress lines, only show the last one (for preprocessing only)
        if task_type == "preprocessing":
            # Use groupby with custom key function to group the download progress lines in order
            def key_func(line):
                return "Downloading: " in line

            log_lines = [
                list(group)[-1] for _, group in groupby(log_lines, key=key_func)
            ]

        # Write filtered log to file, preserving new log output since last read
        tmp_log_path = os.path.join(Path(log_path).parent, ".tmp-log")
        cmd = (f"tail -c {len(full_log_content)} {log_path} > {tmp_log_path} && "
               f"cat {tmp_log_path} > {log_path} && rm {tmp_log_path}")

        full_log_content = "\n".join(log_lines)
        with open(tmp_log_path, "w") as log_file:
            log_file.write(full_log_content)

        subprocess.run(cmd, shell=True, check=True)

        return full_log_content

    def get_log_update(self, last_log_offset=0, full_log=False) -> Tuple[str, int]:
        """Returns the contents of the log file for a task since the last update"""
        log_file_path = os.path.join(self._task_dir, ".log")
        log_content = ""
        log_read_success = False
        if os.path.isfile(log_file_path):
            try:
                full_log_content = TaskRunner._filter_log(log_file_path, self.task_type)
                log_read_success = True
            except UnicodeDecodeError:
                pass

        if log_read_success:
            self._full_log_content = full_log_content
        else:
            full_log_content = self._full_log_content
            with open(log_file_path, "w") as log_file:
                log_file.write(full_log_content)

        full_log_length = len(full_log_content)
        log_content_max_len = 50000
        if (full_log_length == last_log_offset) and not full_log:
            # No new log content since last update
            log_content=""
        elif (full_log_length > last_log_offset) and not full_log:
            # Return just the new part of the log since the last request
            offset = max(0, full_log_length - log_content_max_len)
            offset = max(last_log_offset, offset)
            log_content = full_log_content[offset:]
        else:
            # Return full log content
            log_content = full_log_content[-log_content_max_len:]

        last_log_offset = full_log_length

        return log_content, last_log_offset

    def add_input_file(self, file_contents, standard_fname):
        """Saves uploaded file for a task and returns the save path."""

        if standard_fname not in self._input_filenames:
            raise InvalidTaskInputFile(
                "Failed to add input file: "
                f"Filename '{standard_fname}' is invalid for task type: {self.task_type}. "
                f"Input file must be one of {self._input_filenames}")

        save_path = os.path.join(self._task_dir, "input", standard_fname)
        if os.path.exists(save_path):
            os.remove(save_path)

        lines = file_contents.split(b'\n')

        print(save_path)
        with open(save_path, "w+", encoding="UTF-8") as user_file:
            print(len(lines))
            for line in lines:
                user_file.write(f"{line.decode('UTF-8')}\n")

        return save_path

    def list_input_files(self) -> Dict[str, str]:
        """
        List base filenames of all input files for a task.
        Returns:
            Dict[str, str]: Dictionary of input filenames and their user-specified names.
        """
        redis_db = Redis()
        input_files = {}
        redis_hash_name = f"{self._task_id}_input_files"
        if redis_db.exists(redis_hash_name):
            input_files = redis_db.hgetall(redis_hash_name)
        return {k.decode(): v.decode() for k, v in input_files.items()}

    def remove_task(self):
        """Remove task directory and associated redis data"""
        redis_db = Redis()
        redis_hash_name = f"{self._task_id}_input_files"
        if redis_db.exists(redis_hash_name):
            all_keyval_pairs = redis_db.hgetall(redis_hash_name)
            all_keys = list(all_keyval_pairs.keys())
            redis_db.hdel(redis_hash_name, *all_keys)
        shutil.rmtree(self._task_dir)
