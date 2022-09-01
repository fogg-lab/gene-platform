# http://flask.pocoo.org/docs/1.0/tutorial/database/
import sqlite3

from flask import Flask, current_app, g
app = Flask(__name__)

def get_db():
    """
    Return a database connection.
    If one does not exist, create one.
    """
    if "db" not in g:
        g.db = sqlite3.connect(
            "sqlite_db", detect_types=sqlite3.PARSE_DECLTYPES
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
    db = get_db()

    with current_app.open_resource("schema.sql") as f:
        db.executescript(f.read().decode("utf8"))

def init_db_command():
    """Clear the existing data and create new tables."""
    with app.app_context():
        init_db()

def init_app(app):
    """Register database functions with the Flask app."""
    app.teardown_appcontext(close_db)
    app.cli.add_command(init_db_command)
