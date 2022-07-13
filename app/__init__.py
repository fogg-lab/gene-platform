# imports for debugging (allow printing to stderr)
#from __future__ import print_function
#import sys

from flask import Flask
from flask_session.__init__ import Session
from app.common.common_blueprint import common_bp
from app.analysis.analysis_blueprint import analysis_bp
from app.batch_correction.batch_correction_blueprint import batch_correction_bp
from app.rnaseq_correlation.correlation_blueprint import correlation_bp

def init_app():
    """Initialize the core application."""
    app = Flask(__name__, instance_relative_config=False)
    app.config.from_pyfile("config.py")

    sess = Session()
    sess.init_app(app)

    with app.app_context():

        # Register Blueprints
        app.register_blueprint(common_bp)
        app.register_blueprint(analysis_bp)
        app.register_blueprint(batch_correction_bp)
        app.register_blueprint(correlation_bp)

        return app
