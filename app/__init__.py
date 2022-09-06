import sqlite3
from flask import Flask
from flask_session.__init__ import Session
from flask_login import LoginManager
from app.db.db import init_db_command
from app.models.user import User
from app.blueprints.common import common_bp
from app.blueprints.analysis import analysis_bp
from app.blueprints.batch_correction import batch_correction_bp
from app.blueprints.auth import auth_bp
from app.blueprints.correlation import correlation_bp
from app.blueprints.normalization import normalization_bp
from app.blueprints.preprocessing import preprocessing_bp

def init_app():
    """Initialize the app."""
    app = Flask(__name__, instance_relative_config=False)
    app.config.from_pyfile("config.py")

    sess = Session()
    sess.init_app(app)

    login_manager = LoginManager()
    login_manager.init_app(app)
    @login_manager.user_loader
    def load_user(user):
        return User.get(user)

    try:
        init_db_command()
    except sqlite3.OperationalError:
        # Assume it has already been created
        pass

    with app.app_context():
        # Register Blueprints
        app.register_blueprint(main_bp)
        app.register_blueprint(auth_bp)
        app.register_blueprint(analysis_bp)
        app.register_blueprint(batch_correction_bp)
        app.register_blueprint(correlation_bp)
        app.register_blueprint(normalization_bp)
        app.register_blueprint(preprocessing_bp)

        return app
