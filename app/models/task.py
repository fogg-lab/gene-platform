from __future__ import annotations
import os
import random
import string
from typing import List
from zipfile import ZipFile
from redis import Redis
from rq import Queue
from flask import current_app
from flask_login import current_user

from app.db.db import get_db
from app.task_utils.analysis import AnalysisRunner
from app.task_utils.batch_correction import BatchCorrectionRunner
from app.task_utils.preprocessing import PreprocessingRunner
from app.task_utils.correlation import CorrelationRunner
from app.task_utils.normalization import NormalizationRunner
from app.exceptions import TaskNotFound, InvalidTaskType


class Task:
    """Models class for preparing and queueing user tasks to run in the background."""

    def __init__(self, task_id, user_id, task_type, status, created_at, updated_at):
        self.task_id = task_id
        self.user_id = user_id
        self.task_type = task_type
        self.status = status
        self.created_at = created_at
        self.updated_at = updated_at

    @staticmethod
    def get(task_id: str) -> Task:
        """Retrieve user task by id from database"""
        task = None
        if not isinstance(task_id, str):
            raise TypeError("task_id must be a string")
        if len(task_id) == 16:
            db = get_db()
            user_id = current_user.id
            task = db.execute(
                "SELECT * FROM task WHERE id = ? AND user_id = ?", (task_id, user_id)
            ).fetchone()
        return None if task is None else Task(*task)

    @staticmethod
    def get_log_update(task_id, last_log_offset=0, full_log=False):
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
        log_update = task_runner.get_log_update(last_log_offset, full_log)

        return log_update

    @staticmethod
    def _get_runner(task_id):
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
    def _get_dir(task_id):
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
        return task_dir

    @staticmethod
    def _get_type(task_id):
        """Get task type by id"""
        task = Task.get(task_id)
        return task.task_type

    @staticmethod
    def get_user_tasks():
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
    def submit(task_id):
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
            # queue task
            q = Queue(connection=Redis())
            q.enqueue(
                Task.__execute_task,
                args = (task_id, task_type)
            )
        except TaskNotFound as exc:
            raise TaskNotFound("User task not found, so it could not be submitted.") from exc

    @staticmethod
    def create(task_type, status="pending"):
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
        db.execute(
            "INSERT INTO task (id, user_id, task_type, status) VALUES (?, ?, ?, ?)",
            (new_task_id, user_id, task_type, status)
        )
        db.commit()
        return new_task_id

    @staticmethod
    def delete(task_id):
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
    def delete_input_file(task_id, filename):
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
    def list_input_files(task_id):
        """
        List base filenames of all input files for a task.
        Args:
            task_id (str): The task id.
        Returns:
            list[str]: List of base filenames.
        """
        runner = Task._get_runner(task_id)
        input_file_basenames = runner.list_input_files()
        return input_file_basenames

    @staticmethod
    def get_config(task_id):
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
    def add_input_file(task_id, file_contents, standard_filename, user_filename):
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
        save_path = runner.add_input_file(file_contents, standard_filename)
        # Perform input validation
        status = runner.update_task()
        if "errors" in status:
            os.remove(save_path)
        else:
            redis_db = Redis()
            redis_db.hset(f"{task_id}_input_files", standard_filename, user_filename)
        return status

    @staticmethod
    def configure(task_id, config) -> dict:
        """
        Save configuration for a task.
        Args:
            task_id (str): The task id.
            config (dict): The configuration to save.
        Returns:
            dict: status
        """
        runner = Task._get_runner(task_id)
        status = runner.validate_config(config)
        if len(status.get("errors", [])) == 0:
            runner.save_config(config)
        return status

    @staticmethod
    def validate_task(task_id):
        """
        Validate task input files, including the config file.
        Args:
            task_id (str): The task id.
        Returns:
            dict: status
        """
        runner = Task._get_runner(task_id)
        status = runner.validate_task()
        return status

    @staticmethod
    def __execute_task(task_id):
        """Begin execution of a task."""
        # Get task runner
        runner = Task._get_runner
        # Notify task started
        Task.__notify_task_started(task_id)
        # Run task
        status_msg = runner.execute_task()
        # Notify task completed
        Task.__notify_task_completed(task_id, status_msg)

    @staticmethod
    def __notify_task_started(task_id):
        """Update task status to 'Started'"""
        db = get_db()
        db.execute(
            "UPDATE task SET status = 'Started' WHERE id = ?",
            (task_id,)
        )
        db.commit()

    @staticmethod
    def __notify_task_completed(task_id, status_msg):
        """Update task status to 'Completed'"""
        db = get_db()
        if len(status_msg) > 0:
            status_msg = f"Completed: {status_msg}"
        else:
            status_msg = "Completed"
        db.execute(
            f"UPDATE task SET status = '{status_msg}' WHERE task_id = ?",
            (task_id,)
        )
        db.commit()

    @staticmethod
    def create_task_zip(task_id):
        """Creates a zip file of task output given a task id"""

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

        zip_name = os.path.join(output_dir, task_result_names.get(task_type, "results"))

        with ZipFile(zip_name, "w") as zip_file:
            for file in os.listdir(output_dir):
                zip_file.write(file)

        return zip_name

    @staticmethod
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
