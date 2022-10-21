import os
import random
import string
from redis import Redis
from rq import Queue
from flask import current_app
from flask_login import current_user
from app.db.db import get_db
from app.job_utils.analysis import AnalysisRunner
from app.job_utils.batch_correction import BatchCorrectionRunner
from app.job_utils.correlation import CorrelationRunner
from app.job_utils.normalization import NormalizationRunner
from app.job_utils.preprocessing import PreprocessingRunner
from app.exceptions import JobNotFound

class Job:
    """Models class for preparing and queueing user jobs to run in the background."""

    def __init__(self, job_id, job_type, user_id, status, created_at, updated_at):
        self.job_id = job_id
        self.user_id = user_id
        self.status = status
        self.created_at = created_at
        self.updated_at = updated_at
        self.job_type = job_type

    @staticmethod
    def get(job_id):
        """Retrieve user job by id from database"""
        job = None
        if job_id is not None and len(job_id == 16):
            db = get_db()
            user_id = current_user.user_id
            job = db.execute(
                "SELECT * FROM job WHERE id = ? AND user_id = ?", (job_id, user_id)
            ).fetchone()
        if job is None:
            return None
        return Job(
            job_id=job[0],
            user_id=job[1],
            job_type=job[2],
            status=job[3],
            created_at=job[4],
            updated_at=job[5])

    @staticmethod
    def _get_runner(job_id):
        """Retrieve runner for a job by id"""
        job_dir = Job._get_dir(job_id)
        job_type = Job._get_type(job_id)
        if job_type == "analysis":
            return AnalysisRunner(job_id, job_dir)
        elif job_type == "batch_correction":
            return BatchCorrectionRunner(job_id, job_dir)
        elif job_type == "correlation":
            return CorrelationRunner(job_id, job_dir)
        elif job_type == "normalization":
            return NormalizationRunner(job_id, job_dir)
        elif job_type == "preprocessing":
            return PreprocessingRunner(job_id, job_dir)

    @staticmethod
    def _get_dir(job_id):
        """
        Return absolute path to the directory for a job associated with a user.
        Create the directory if it doesn't exist.
        """
        try:
            Job.get(job_id)
        except JobNotFound as exc:
            raise JobNotFound(
                "User job not found, so no directory was returned.") from exc
        user_jobs_path = current_app.config["USER_JOBS_PATH"]
        job_dir = os.path.join(user_jobs_path, job_id)
        os.makedirs(job_dir, exist_ok=True)
        return job_dir

    @staticmethod
    def _get_type(job_id):
        """Get job type by id"""
        job = Job.get(job_id)
        return job.job_type

    @staticmethod
    def get_user_jobs():
        """Retrieve all jobs for the current user"""
        user_id = current_user.user_id
        db = get_db()
        jobs = db.execute(
            "SELECT * FROM job WHERE user_id = ?", (user_id,)
        ).fetchall()
        jobs = [Job(
            job_id=job[0],
            user_id=job[1],
            job_type=job[2],
            status=job[3],
            created_at=job[4],
            updated_at=job[5]
        ) for job in jobs]
        return jobs

    @staticmethod
    def submit(job_id):
        """Update status of user job to 'ready' and queue it for processing."""
        try:
            Job.get(job_id)
            db = get_db()
            db.execute(
                "UPDATE job SET status = 'ready' WHERE id = ?",
                (job_id,)
            )
            # get job type from db
            job_type = db.execute(
                "SELECT job_type FROM job WHERE id = ?",
                (job_id,)
            ).fetchone()[0]
            # queue job
            q = Queue(connection=Redis())
            q.enqueue(
                Job.__start_job,
                args = (job_id, job_type)
            )
        except JobNotFound as exc:
            raise JobNotFound("User job not found, so it could not be submitted.") from exc

    @staticmethod
    def create(job_type, status="pending"):
        """
        Add a job to the database.
        Args:
            job_type (string): The type of job to create.
            status (string): The initial status of job (default: pending).
        """
        db = get_db()
        new_job_id = ''.join(random.choice(string.ascii_lowercase + string.digits)
                             for _ in range(16))
        user_id = current_user.user_id
        db.execute(
            "INSERT INTO job (id, user_id, job_type, status) VALUES (?, ?, ?, ?)",
            (new_job_id, user_id, job_type, status)
        )
        db.commit()
        return new_job_id

    @staticmethod
    def delete(job_id):
        """Delete job from db"""
        job_found = False
        try:
            Job.get(job_id)
            job_found = True
        except JobNotFound as exc:
            raise JobNotFound(
                "User job not found, so it could not be deleted.") from exc
        if job_found:
            db = get_db()
            db.execute(
                "DELETE FROM job WHERE id = ?", (job_id,)
            )
            db.commit()

    @staticmethod
    def delete_input_file(job_id, filename):
        """
        Remove an input file for a job.
        Args:
            job_id (str): The job id.
            filename (str): The name of the input file to be deleted.
        """
        job_dir = Job._get_dir(job_id)
        input_file_path = os.path.join(job_dir, filename)
        if os.path.isfile(input_file_path):
            os.remove(input_file_path)

    @staticmethod
    def list_input_files(job_id):
        """
        List base filenames of all input files for a job.
        Args:
            job_id (str): The job id.
        Returns:
            list[str]: List of base filenames.
        """
        runner = Job._get_runner(job_id)
        input_file_basenames = runner.list_input_files()
        return input_file_basenames

    @staticmethod
    def add_input_file(job_id, file_contents, standard_filename, user_filename):
        """
        Saves uploaded file for a job and perform input validation
        Args:
            file_contents (bytes): File contents as a bytes stream.
            standard_filename (str): Filename expected by the job runner.
            user_filename (str): Original name of the user-uploaded file.
        Returns:
            dict: status
        """
        runner = Job._get_runner(job_id)
        save_path = runner.add_input_file(file_contents, standard_filename)
        # Perform input validation
        status = runner.update_job()
        if "errors" in status:
            os.remove(save_path)
        else:
            redis_db = Redis()
            redis_db.hset(f"{job_id}_input_files", standard_filename, user_filename)

    @staticmethod
    def configure(job_id, config):
        """
        Save configuration for a job.
        Args:
            job_id (str): The job id.
            config (dict): The configuration to save.
        """
        runner = Job._get_runner(job_id)
        status = runner.validate_config(config)
        if len(status["errors"]) == 0:
            runner.save_config(config)

    @staticmethod
    def __start_job(job_id):
        """Begin execution of a job."""
        # Get job runner
        runner = Job._get_runner
        # Notify job started
        Job.__notify_job_started(job_id)
        # Run job
        status_msg = runner.start_job()
        # Notify job completed
        Job.__notify_job_completed(job_id, status_msg)

    @staticmethod
    def __notify_job_started(job_id):
        """Update job status to 'Started'"""
        db = get_db()
        db.execute(
            "UPDATE job SET status = 'Started' WHERE id = ?",
            (job_id,)
        )
        db.commit()

    @staticmethod
    def __notify_job_completed(job_id, status_msg):
        """Update job status to 'Completed'"""
        db = get_db()
        if len(status_msg) > 0:
            status_msg = f"Completed: {status_msg}"
        else:
            status_msg = "Completed"
        db.execute(
            f"UPDATE job SET status = '{status_msg}' WHERE job_id = ?",
            (job_id,)
        )
        db.commit()
