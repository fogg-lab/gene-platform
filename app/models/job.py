import os
import random
import string
from redis import Redis
from rq import Queue
from flask import current_app
from flask_login import current_user
from app.db.db import get_db
from app.job_runner.job_runner import run_job
from app.exceptions import JobNotFoundError

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
            raise JobNotFoundError
        return Job(
            job_id=job[0],
            user_id=job[1],
            job_type=job[2],
            status=job[3],
            created_at=job[4],
            updated_at=job[5])

    @staticmethod
    def get_dir(job_id):
        """
        Return absolute path to the directory for a job associated with a user.
        Create the directory if it doesn't exist.
        """
        try:
            Job.get(job_id)
        except JobNotFoundError as exc:
            raise JobNotFoundError(
                "User job not found, so no directory was returned.") from exc
        user_jobs_path = current_app.config["USER_JOBS_PATH"]
        job_dir = os.path.join(user_jobs_path, job_id)
        os.makedirs(job_dir, exist_ok=True)
        return job_dir

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
    def submit_job(job_id):
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
        except JobNotFoundError as exc:
            raise JobNotFoundError(
                "User job not found, so it could not be submitted.") from exc

    @staticmethod
    def create(job_type, status="pending"):
        """Add job to db"""
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
        except JobNotFoundError as exc:
            raise JobNotFoundError(
                "User job not found, so it could not be deleted.") from exc
        if job_found:
            db = get_db()
            db.execute(
                "DELETE FROM job WHERE id = ?", (job_id,)
            )
            db.commit()

    @staticmethod
    def __start_job(job_id, job_type):
        """Start job"""
        # Get job directory
        job_dir = Job.get_dir(job_id)
        # Notify job started
        Job.__notify_job_started(job_id)
        # Run job
        status_msg = run_job(job_dir, job_type)
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
