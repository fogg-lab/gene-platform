from __future__ import annotations
import os
import random
import string
from pathlib import Path
from typing import List, Tuple
from zipfile import ZipFile, ZIP_DEFLATED
from glob import glob
from functools import wraps
from redis import Redis
from flask import current_app
from flask_login import current_user

from app.db.db import get_db
from app.task_utils.task_runner import TaskRunner
from app.task_utils.analysis import AnalysisRunner
from app.task_utils.batch_correction import BatchCorrectionRunner
from app.task_utils.preprocessing import PreprocessingRunner
from app.task_utils.correlation import CorrelationRunner
from app.task_utils.normalization import NormalizationRunner
from app.exceptions import TaskNotFound, InvalidTaskType
from app.helper import StatusDict


def require_task_id_correct_format(task_method):
    """Decorator that raises an exception if task_id is not a 16-char long alphanumeric string."""
    @wraps(task_method)

    def check_task_id(*args, **kwargs):
        if "task_id" in kwargs:
            task_id = kwargs["task_id"]
        elif len(args) == 0:
            raise ValueError("task_id not found in function arguments")
        else:
            task_id = args[0]
        if not isinstance(task_id, str):
            raise TypeError("task_id must be a string")
        if len(task_id) != 16:
            raise ValueError("task_id must be 16 characters")
        if not task_id.isalnum():
            raise ValueError("task_id must be alphanumeric")

        return task_method(*args, **kwargs)

    return check_task_id

class Task:
    """Models class for preparing and queueing user tasks to run in the background."""

    def __init__(self, task_id, user_id, task_type, status, created_at, updated_at):
        if not task_id.isalnum() or len(task_id) != 16:
            raise ValueError("task_id must be 16 characters and alphanumeric")
        self.task_id = task_id
        self.user_id = user_id
        self.task_type = task_type
        self.status = status
        self.created_at = created_at
        self.updated_at = updated_at

    @staticmethod
    @require_task_id_correct_format
    def get(task_id: str) -> Task:
        """Retrieve user task by id from database"""
        task = None
        if len(task_id) == 16:
            db = get_db()
            user_id = current_user.id
            cur = db.cursor()
            task = cur.execute(
                "SELECT * FROM task WHERE id = ? AND user_id = ?", (task_id, user_id)
            ).fetchone()
        return None if task is None else Task(*task)

    @staticmethod
    @require_task_id_correct_format
    def get_log_update(task_id, last_log_offset=0, full_log=False) -> Tuple[str, int]:
        """
        Returns the contents of the log file for a task
        Args:
            task_id (str): 16-character alphanumeric task id
            last_log_offset (int): Length of full log at the last log update. Defaults to 0.
            full_log (bool): Whether to return the full log (instead of offset). Defaults to False.
        """

        if last_log_offset == 0:
            full_log = True
        task_runner = Task._get_runner(task_id)
        log_content, last_log_offset = task_runner.get_log_update(last_log_offset, full_log)

        return log_content, last_log_offset

    @staticmethod
    @require_task_id_correct_format
    def _get_runner(task_id) -> TaskRunner:
        """Retrieve runner for a task by id"""
        task_dir = Task._get_dir(task_id)
        task_type = Task._get_type(task_id)
        if task_type == "analysis":
            return AnalysisRunner(task_id, task_dir)
        elif task_type == "batch_correction":
            return BatchCorrectionRunner(task_id, task_dir)
        elif task_type == "correlation":
            return CorrelationRunner(task_id, task_dir)
        elif task_type == "normalization":
            return NormalizationRunner(task_id, task_dir)
        elif task_type == "preprocessing":
            return PreprocessingRunner(task_id, task_dir)
        else:
            raise InvalidTaskType(f"Unrecognized task type: {task_type}")

    @staticmethod
    @require_task_id_correct_format
    def _get_dir(task_id: str) -> str:
        """
        Return absolute path to the directory for a task associated with a user.
        Create the directory if it doesn't exist.
        """
        try:
            Task.get(task_id)
        except TaskNotFound as exc:
            raise TaskNotFound(
                "User task not found, so no directory was returned.") from exc
        user_tasks_path = current_app.config["USER_TASKS_PATH"]
        task_dir = os.path.join(user_tasks_path, task_id)
        os.makedirs(task_dir, exist_ok=True)
        input_dir = os.path.join(task_dir, "input")
        output_dir = os.path.join(task_dir, "output")
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        return task_dir

    @staticmethod
    @require_task_id_correct_format
    def _get_type(task_id) -> str:
        """Get task type by id"""
        task = Task.get(task_id)
        return task.task_type

    @staticmethod
    def get_user_tasks() -> List[Task]:
        """Retrieve all tasks for the current user"""
        user_id = current_user.id
        db = get_db()
        tasks = db.execute(
            "SELECT * FROM task WHERE user_id = ?", (user_id,)
        ).fetchall()
        tasks = [Task(
            task_id=task[0],
            user_id=task[1],
            task_type=task[2],
            status=task[3],
            created_at=task[4],
            updated_at=task[5]
        ) for task in tasks]
        return tasks

    @staticmethod
    @require_task_id_correct_format
    def submit(task_id) -> str:
        """Update status of user task to 'ready' and queue it for processing."""
        try:
            Task.get(task_id)
            db = get_db()
            db.execute(
                "UPDATE task SET status = 'ready' WHERE id = ?",
                (task_id,)
            )
            # get task type from db
            task_type = db.execute(
                "SELECT task_type FROM task WHERE id = ?",
                (task_id,)
            ).fetchone()[0]
            runner = Task._get_runner(task_id)
            Task._update_task_status(task_id, "queued")
            runner.execute_task()
        except TaskNotFound as exc:
            raise TaskNotFound("User task not found, so it could not be submitted.") from exc

        return ""

    @staticmethod
    def create(task_type, status="pending") -> Task:
        """
        Add a task to the database.
        Args:
            task_type (string): The type of task to create.
            status (string): The initial status of task (default: pending).
        """
        db = get_db()
        new_task_id = ''.join(random.choice(string.ascii_lowercase + string.digits)
                             for _ in range(16))
        user_id = current_user.id
        cur = db.cursor()
        cur.execute(
            "INSERT INTO task (id, user_id, task_type, status) VALUES (?, ?, ?, ?)",
            (new_task_id, user_id, task_type, status)
        )
        db.commit()
        return new_task_id

    @staticmethod
    @require_task_id_correct_format
    def delete(task_id) -> None:
        """Delete task from db and remove task directory"""
        try:
            runner = Task._get_runner(task_id)
        except TaskNotFound as exc:
            raise TaskNotFound("User task not found, so it could not be deleted.") from exc

        runner.remove_task()

        db = get_db()
        db.execute(
            "DELETE FROM task WHERE id = ?", (task_id,)
        )
        db.commit()

    @staticmethod
    @require_task_id_correct_format
    def delete_input_file(task_id, filename) -> None:
        """
        Remove an input file for a task.
        Args:
            task_id (str): The task id.
            filename (str): The name of the input file to be deleted.
        """
        task_dir = Task._get_dir(task_id)
        input_file_path = os.path.join(task_dir, filename)
        if os.path.isfile(input_file_path):
            os.remove(input_file_path)

    @staticmethod
    @require_task_id_correct_format
    def list_input_files(task_id) -> List[str]:
        """
        List base filenames of all input files for a task.
        Args:
            task_id (str): The task id.
        Returns:
            List[str]: List of base filenames.
        """
        runner = Task._get_runner(task_id)
        input_file_basenames = runner.list_input_files()
        return input_file_basenames

    @staticmethod
    @require_task_id_correct_format
    def get_config(task_id) -> dict:
        """
        Get task configuration.
        Args:
            task_id (str): The task id.
        Returns:
            dict: The configuration parameters for the task.
        """
        runner = Task._get_runner(task_id)
        config = runner.get_config()
        return config

    @staticmethod
    @require_task_id_correct_format
    def add_input_file(task_id, file_contents, standard_filename, user_filename) -> StatusDict:
        """
        Saves uploaded file for a task and perform input validation
        Args:
            file_contents (bytes): File contents as a bytes stream.
            standard_filename (str): Filename expected by the task runner.
            user_filename (str): Original name of the user-uploaded file.
        Returns:
            dict: status
        """
        runner = Task._get_runner(task_id)
        #save_path = runner.add_input_file(file_contents, standard_filename)
        runner.add_input_file(file_contents, standard_filename)
        # Perform input validation
        #status = runner.update_task()
        status=StatusDict()
        #if "errors" in status:
            #os.remove(save_path)
        #else:
        redis_db = Redis()
        redis_db.hset(f"{task_id}_input_files", standard_filename, user_filename)
        return status

    @staticmethod
    @require_task_id_correct_format
    def configure(task_id, config) -> StatusDict:
        """
        Save configuration for a task.
        Args:
            task_id (str): The task id.
            config (dict): The configuration to save.
        Returns:
            StatusDict: A status message and errors, if any.
        """
        runner = Task._get_runner(task_id)
        status = runner.validate_config(config)
        if len(status.get("errors", [])) == 0:
            runner.save_config(config)
        return status

    @staticmethod
    @require_task_id_correct_format
    def validate_task(task_id) -> StatusDict:
        """
        Validate task input files, including the config file.
        Args:
            task_id (str): The task id.
        Returns:
            StatusDict: A status message and errors, if any.
        """
        runner = Task._get_runner(task_id)
        status = runner.validate_task()
        return status

    @staticmethod
    @require_task_id_correct_format
    def _update_task_status(task_id, status) -> None:
        """Update task status to 'Started'"""
        db = get_db()
        db.execute("UPDATE task SET status = ? WHERE id = ?", (status, task_id))
        db.commit()

    @staticmethod
    @require_task_id_correct_format
    def get_status(task_id) -> str:
        """Get status of a task."""
        db = get_db()
        status = db.execute(
            "SELECT status FROM task WHERE id = ?",
            (task_id,)
        ).fetchone()[0]
        status = "" if status is None else status
        return status

    @staticmethod
    @require_task_id_correct_format
    def create_task_zip(task_id) -> str:
        """Creates a zip file of task output given a task id and returns the filepath."""

        task_dir = Task._get_dir(task_id)
        task_type = Task._get_type(task_id)

        output_dir = os.path.join(task_dir, "output")

        task_result_names = {
            "correlation": "correlation_plots",
            "batch_correction": "batch_corrected_counts",
            "normalization": "normalized_counts",
            "preprocessing": "preprocessed_data",
            "analysis": "analysis_results"
        }

        zip_path = os.path.join(output_dir, task_result_names.get(task_type, "results") + ".zip")
        if os.path.isfile(zip_path):
            os.remove(zip_path)

        out_filepaths = [fp for fp in glob(os.path.join(output_dir, "*")) if os.path.isfile(fp)]

        if len(out_filepaths) == 0:
            return ""

        with ZipFile(zip_path, mode="w", compression=ZIP_DEFLATED) as zip_file:
            for filepath in out_filepaths:
                zip_file.write(filepath, arcname=Path(filepath).name)

        return zip_path

    @staticmethod
    @require_task_id_correct_format
    def get_output_filepath(task_id, standard_filename) -> str:
        """
        Get the path to a task output file.
        Args:
            task_id (str): The task id.
            standard_filename (str): The filename expected by the task runner.
        Returns:
            str: The path to the output file.
        """
        runner = Task._get_runner(task_id)
        output_filepath = runner.get_output_filepath(standard_filename)
        return output_filepath

    @staticmethod
    @require_task_id_correct_format
    def get_all_output_filepaths(task_id) -> List[str]:
        """
        Get the paths to all task output files.
        Args:
            task_id (str): The task id.
        Returns:
            list[str]: The paths to the output files.
        """
        runner = Task._get_runner(task_id)
        output_filepaths = runner.get_all_output_filepaths()
        return output_filepaths
