"""Standalone task execution function for rq to run asynchonously, outside of the app context."""
import subprocess
import sqlite3
from pathlib import Path
import os.path


def execute_job_async(task_id, log_path, task_type, cmd):
    """Run a task, then update the task status in the userdata database."""
    expected_db_filepath = Path(os.path.abspath(__file__)).parent.parent / "db" / "userdata.db"
    if os.path.isfile(expected_db_filepath):
        db = sqlite3.connect(
            expected_db_filepath, detect_types=sqlite3.PARSE_DECLTYPES
        )
        db.row_factory = sqlite3.Row
        cursor = db.cursor()
        cursor.execute("UPDATE task SET status = 'running' WHERE id = ?", (task_id,))
        db.commit()
        db.close()
    else:
        raise FileNotFoundError(f"Database file not found: {expected_db_filepath}")
    with open(log_path, "w") as log_file:
        log_file.write(f"Starting {task_type}...\n")
    try:
        subprocess.check_call(cmd+" >>{log_path} 2>&1", shell=True)
        if os.path.isfile(expected_db_filepath):
            db = sqlite3.connect(
                expected_db_filepath, detect_types=sqlite3.PARSE_DECLTYPES
            )
            db.row_factory = sqlite3.Row
            cursor = db.cursor()
            cursor.execute("UPDATE task SET status = 'completed' WHERE id = ?", (task_id,))
            db.commit()
            db.close()
        else:
            raise FileNotFoundError(f"Database file not found: {expected_db_filepath}")
    except subprocess.CalledProcessError as err:
        if os.path.isfile(expected_db_filepath):
            db = sqlite3.connect(
                expected_db_filepath, detect_types=sqlite3.PARSE_DECLTYPES
            )
            db.row_factory = sqlite3.Row
            cursor = db.cursor()
            cursor.execute("UPDATE task SET status = 'failed' WHERE id = ?", (task_id,))
            db.commit()
            db.close()
        else:
            raise FileNotFoundError(f"Database file not found: {expected_db_filepath}") from err
        raise err
