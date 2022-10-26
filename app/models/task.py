import os
import random
import string
from redis import Redis
from rq import Queue
from flask import current_app
from flask_login import current_user
from app.db.db import get_db
from app.task_utils.analysis import AnalysisRunner
from app.task_utils.batch_correction import BatchCorrectionRunner
from app.task_utils.correlation import CorrelationRunner
from app.task_utils.normalization import NormalizationRunner
from app.task_utils.preprocessing import PreprocessingRunner
from app.exceptions import TaskNotFound

class Task:
    """Models class for preparing and queueing user tasks to run in the background."""

    def __init__(self, task_id, task_type, user_id, status, created_at, updated_at):
        self.task_id = task_id
        self.user_id = user_id
        self.status = status
        self.created_at = created_at
        self.updated_at = updated_at
        self.task_type = task_type

    @staticmethod
    def get(task_id):
        """Retrieve user task by id from database"""
        task = None
        if task_id is not None and len(task_id == 16):
            db = get_db()
            user_id = current_user.user_id
            task = db.execute(
                "SELECT * FROM task WHERE id = ? AND user_id = ?", (task_id, user_id)
            ).fetchone()
        if task is None:
            return None
        return Task(
            task_id=task[0],
            user_id=task[1],
            task_type=task[2],
            status=task[3],
            created_at=task[4],
            updated_at=task[5])

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
        user_id = current_user.user_id
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
                Task.__start_task,
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
        user_id = current_user.user_id
        db.execute(
            "INSERT INTO task (id, user_id, task_type, status) VALUES (?, ?, ?, ?)",
            (new_task_id, user_id, task_type, status)
        )
        db.commit()
        return new_task_id

    @staticmethod
    def delete(task_id):
        """Delete task from db"""
        task_found = False
        try:
            Task.get(task_id)
            task_found = True
        except TaskNotFound as exc:
            raise TaskNotFound(
                "User task not found, so it could not be deleted.") from exc
        if task_found:
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

    @staticmethod
    def configure(task_id, config):
        """
        Save configuration for a task.
        Args:
            task_id (str): The task id.
            config (dict): The configuration to save.
        """
        runner = Task._get_runner(task_id)
        status = runner.validate_config(config)
        if len(status["errors"]) == 0:
            runner.save_config(config)

    @staticmethod
    def __start_task(task_id):
        """Begin execution of a task."""
        # Get task runner
        runner = Task._get_runner
        # Notify task started
        Task.__notify_task_started(task_id)
        # Run task
        status_msg = runner.start_task()
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
