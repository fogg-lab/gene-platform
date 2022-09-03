import os
import random
import string
from redis import Redis
from rq import Queue
from flask import current_app
from app.db.db import get_db
from app.job_tools import run_job, prepare_job

class Job:
    """Models class for user jobs"""

    accepted_job_types = ["batch_correction", "analysis", "preprocessing", "correlation"]

    def __init__(self, job_id, job_type, user_id, status, created_at, updated_at):
        self.job_id = job_id
        self.user_id = user_id
        self.status = status
        self.created_at = created_at
        self.updated_at = updated_at
        if job_type not in Job._job_modules:
            raise ValueError(f"Invalid job type: {job_type}.\n"
                             f"Job type be one of {list(Job._job_modules.keys())}")
        self.job_type = job_type

    @staticmethod
    def get_job_path(job_id):
        """
        Return absolute path to the job directory.
        Create the directory if it doesn't exist
        """
        user_jobs_path = current_app.config["USER_JOBS_PATH"]
        job_path = os.path.join(user_jobs_path, job_id)
        os.makedirs(job_path, exist_ok=True)
        return job_path

    @staticmethod
    def queue_job(job_id):
        """Update job status to 'ready' and queue the job"""
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

    @staticmethod
    def __notify_job_started(job_id):
        """Update job status to 'started'"""
        db = get_db()
        db.execute(
            "UPDATE job SET status = 'started' WHERE id = ?",
            (job_id,)
        )
        db.commit()


    @staticmethod
    def __notify_job_completed(job_id):
        """Update job status to 'completed'"""
        db = get_db()
        db.execute(
            "UPDATE job SET status = 'completed' WHERE job_id = ?",
            (job_id,)
        )
        db.commit()

    @staticmethod
    def __start_job(job_id, job_type):
        """Start job"""
        # Get job directory
        job_dir = Job.get_job_path(job_id)
        # Notify job started
        Job.__notify_job_started(job_id)
        # Run job script
        Job._job_modules[job_type].start_job(job_dir)
        # Notify job completed
        Job.__notify_job_completed(job_id)

    @staticmethod
    def get(job_id):
        """Retrieve job by id from database"""
        db = get_db()
        job = db.execute(
            "SELECT * FROM job WHERE id = ?", (job_id,)
        ).fetchone()
        if not job:
            return None
        job = Job(
            job_id=job[0],
            user_id=job[1],
            job_type=job[2],
            status=job[3],
            created_at=job[4],
            updated_at=job[5]
        )
        return job

    @staticmethod
    def get_user_jobs(user_id):
        """Retrieve all jobs for a user"""
        db = get_db()
        jobs = db.execute(
            "SELECT * FROM job WHERE user_id = ?", (user_id,)
        ).fetchall()
        if not jobs:
            return None
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
    def create(status="pending"):
        """Add job to db"""
        db = get_db()
        new_job_id = ''.join(random.choice(string.ascii_lowercase + string.digits)
                             for _ in range(16))
        db.execute(
            "INSERT INTO job (id, status) VALUES (?, ?)",
            (new_job_id, status)
        )
        db.commit()
        return new_job_id

    @staticmethod
    def delete(job_id):
        """Delete job from db"""
        db = get_db()
        db.execute(
            "DELETE FROM job WHERE id = ?", (job_id,)
        )
        db.commit()
