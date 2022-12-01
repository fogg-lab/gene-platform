import os
from base64 import b64encode
import subprocess
from datetime import timedelta
from dotenv import load_dotenv

load_dotenv()

TEMPLATES_AUTO_RELOAD = True
REDIS_PORT = os.environ.get("REDIS_PORT", 6379)
TEMPLATES_AUTO_RELOAD = True

### Authentication ###

# ENABLE_GOOGLE_AUTH enables "Sign in with Google" authentication
ENABLE_GOOGLE_AUTH = os.environ.get("ENABLE_GOOGLE_AUTH", False)
GOOGLE_CLIENT_ID = os.environ.get("GOOGLE_CLIENT_ID", None)
GOOGLE_CLIENT_SECRET = os.environ.get("GOOGLE_CLIENT_SECRET", None)
GOOGLE_DISCOVERY_URL = "https://accounts.google.com/.well-known/openid-configuration"
os.environ["GOOGLE_DISCOVERY_URL"] = GOOGLE_DISCOVERY_URL

SECRET_KEY = os.getenv("SESSION_SECRET_KEY", None)
if SECRET_KEY is None:
    os.environ["SESSION_SECRET_KEY"] = b64encode(os.urandom(24)).decode("utf-8")
    SECRET_KEY = os.getenv("SESSION_SECRET_KEY")
    subprocess.Popen([f"dotenv set SESSION_SECRET_KEY {SECRET_KEY} > /dev/null"], shell=True)

SESSION_PERMANENT = True
SESSION_TYPE = "filesystem"
SESSION_COOKIE_NAME = "my_session"
PERMANENT_SESSION_LIFETIME = timedelta(minutes=30)
DEBUG = True

SQLALCHEMY_TRACK_MODIFICATIONS = os.environ.get("SQLALCHEMY_TRACK_MODIFICATIONS", False)
SQLALCHEMY_DATABASE_URI = os.environ.get("SQLALCHEMY_DATABASE_URI", "sqlite:///user.sqlite3")

APP_ROOT_PATH = os.path.dirname(os.path.abspath(__file__))
GDC_DATA_PATH = os.path.join(APP_ROOT_PATH,"data/preprocessing/gdc_data")
GEO_DATA_PATH = os.path.join(APP_ROOT_PATH,"data/preprocessing/geo_data")
PROBEMAP_PATH = os.path.join(GEO_DATA_PATH,"probe_maps")
USER_TASKS_PATH = os.path.join(APP_ROOT_PATH, "data/user_tasks")
SCRIPTS_PATH = os.path.join(APP_ROOT_PATH, "scripts")
