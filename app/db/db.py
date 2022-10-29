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

        user_table_sql, task_table_sql = sql_script.split("/* table separator */")

        g.db = sqlite3.connect(
            db_path, detect_types=sqlite3.PARSE_DECLTYPES
        )
        g.db.row_factory = sqlite3.Row

        cursor = g.db.cursor()

        print(task_table_sql)
        cursor.execute(user_table_sql)
        g.db.commit()
        cursor = g.db.cursor()
        cursor.execute(task_table_sql)

        g.db.close()


def init_db_command():
    """Clear the existing data and create new tables."""
    with app.app_context():
        init_db()

def init_app(app):
    """Register database functions with the Flask app."""
    app.teardown_appcontext(close_db)
    app.cli.add_command(init_db_command)
