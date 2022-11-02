# http://flask.pocoo.org/docs/1.0/tutorial/database/
import os.path
from pathlib import Path
import sqlite3
from flask import Flask, g
app = Flask(__name__)

def get_db():
    """
    Return a database connection.
    If one does not exist, create one.
    """

    expected_db_filepath = Path(os.path.abspath(__file__)).parent / "userdata.db"

    if "db" in g:
        return g.db
    elif os.path.isfile(expected_db_filepath):
        g.db = sqlite3.connect(
            expected_db_filepath, detect_types=sqlite3.PARSE_DECLTYPES
        )
        g.db.row_factory = sqlite3.Row
        return g.db

def close_db(e=None):
    """Close the database connection."""
    db = g.pop("db", None)

    if db is not None:
        db.close()

def init_db():
    """Initialize the database."""

    g.db = get_db()
    if g.db is None:
        db_dir_path = Path(os.path.abspath(__file__)).parent
        db_path = db_dir_path / "userdata.db"
        schema_path = os.path.join(db_dir_path, "schema.sql")

        with open(schema_path, 'r') as sql_file:
            sql_script = sql_file.read()

        g.db = sqlite3.connect(
            db_path, detect_types=sqlite3.PARSE_DECLTYPES
        )
        g.db.row_factory = sqlite3.Row
        for sql_statement in sql_script.split("-- separator --"):
            cursor = g.db.cursor()
            cursor.execute(sql_statement)
            g.db.commit()

        g.db.close()


def init_db_command():
    """Clear the existing data and create new tables."""
    with app.app_context():
        init_db()

def init_app(app):
    """Register database functions with the Flask app."""
    app.teardown_appcontext(close_db)
    app.cli.add_command(init_db_command)

# CREATE TABLE user (
#   id TEXT PRIMARY KEY,
#   name TEXT,
#   email TEXT UNIQUE,
#   is_guest INTEGER DEFAULT 0,
#   created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
# );




# CREATE TABLE task (
#   id TEXT NOT NULL,
#   user_id TEXT NOT NULL,
#   task_type TEXT NOT NULL,
#   status TEXT NOT NULL,
#   created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
#   updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
#   PRIMARY KEY (id),
#   FOREIGN KEY (user_id) REFERENCES user(id)
# );




# CREATE TRIGGER [UpdateLastTime]
#     AFTER UPDATE
#     ON task
#     FOR EACH ROW
#     WHEN NEW.updated_at < OLD.updated_at    --- this avoid infinite loop
# BEGIN
#     UPDATE task SET updated_at=CURRENT_TIMESTAMP WHERE id=OLD.id;
# END;