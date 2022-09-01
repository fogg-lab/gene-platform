import os
from base64 import b64encode
import subprocess
from datetime import timedelta
from dotenv import load_dotenv

load_dotenv()

### Authentication ###

# ENABLE_GOOGLE_AUTH enables "Sign in with Google" authentication
ENABLE_GOOGLE_AUTH = os.environ.get("ENABLE_GOOGLE_AUTH", False)
GOOGLE_CLIENT_ID = os.environ.get("GOOGLE_CLIENT_ID", None)
GOOGLE_CLIENT_SECRET = os.environ.get("GOOGLE_CLIENT_SECRET", None)
GOOGLE_DISCOVERY_URL = "https://accounts.google.com/.well-known/openid-configuration"
os.environ["GOOGLE_DISCOVERY_URL"] = GOOGLE_DISCOVERY_URL

# ENABLE_PROFILE_LOGIN adds the ability for the user to login to a profile with no password
# User selects profile from a list of all profiles, or registers a new profile
ENABLE_PROFILE_LOGIN = os.environ.get("ENABLE_PROFILE_LOGIN", True)

# ENABLE_USERNAME_PASSWORD_LOGIN allows the user to login with a username and password
# when REQUIRE_EMAIL_CONFIRMATION is set, users are required to confirm their email address
# ALLOWED_EMAIL_DOMAINS is a list of email domains that are allowed to register
ENABLE_USERNAME_PASSWORD_LOGIN = os.environ.get("ENABLE_USERNAME_PASSWORD_LOGIN", False)
REQUIRE_EMAIL_CONFIRMATION = os.environ.get("REQUIRE_EMAIL_CONFIRMATION", False)
ALLOWED_EMAIL_DOMAINS = os.environ.get("ALLOWED_EMAIL_DOMAINS", "*").replace(" ", "").split(",")

# ENABLE_GUEST_LOGIN allows the user to continue without logging in
# In this case, the user's session is temporary, uses cookies
# GUEST_SESSION_LIFESPAN specifies how long (in hours) guest user jobs are kept on the server
# Once a job's lifespan is over, the files associated with the job are queued for deletion
ENABLE_GUEST_LOGIN = os.environ.get("ENABLE_GUEST_LOGIN", True)
GUEST_SESSION_LIFESPAN = timedelta(hours=os.environ.get("GUEST_SESSION_LIFESPAN", 3))
GUEST_JOB_LIMIT = os.environ.get("GUEST_JOB_LIMIT", 5)

SECRET_KEY = os.getenv("SESSION_SECRET_KEY", None)
if SECRET_KEY is None:
    os.environ["SESSION_SECRET_KEY"] = b64encode(os.urandom(24)).decode("utf-8")
    SECRET_KEY = os.getenv("SESSION_SECRET_KEY")
    subprocess.Popen([f"dotenv set SESSION_SECRET_KEY {SECRET_KEY} > /dev/null"], shell=True)

SESSION_PERMANENT = True
SESSION_TYPE = "filesystem"
SESSION_COOKIE_NAME = "my_session"
PERMANENT_SESSION_LIFETIME = timedelta(hours=GUEST_SESSION_LIFESPAN)
DEBUG = os.environ.get("DEBUG", False)

# USER_JOB_LIMIT specifies the max number of most recent jobs that are kept for each user
# USER_JOB_LIFESPAN specifies how long (in hours) user jobs are kept on the server
# Once a job's lifespan is over, the files associated with the job are queued for deletion
# Applies to users that sign in with username/password, Google, or with profile login
USER_JOB_LIMIT = os.environ.get("USER_JOB_LIMIT", 10)
USER_JOB_LIFESPAN = timedelta(hours=os.environ.get("USER_JOB_LIFESPAN", 168))

SQLALCHEMY_TRACK_MODIFICATIONS = os.environ.get("SQLALCHEMY_TRACK_MODIFICATIONS", False)
SQLALCHEMY_DATABASE_URI = os.environ.get("SQLALCHEMY_DATABASE_URI", "sqlite:///user.sqlite3")

APP_ROOT_PATH = os.path.dirname(os.path.abspath(__file__))
GDC_DATA_PATH = os.path.join(APP_ROOT_PATH,"data/gdc_data")
USER_JOBS_PATH = os.path.join(APP_ROOT_PATH, "data/user_jobs")
RSCRIPTS_PATH = os.path.join(APP_ROOT_PATH, "rscripts")
JOB_SCRIPTS_PATH = os.path.join(APP_ROOT_PATH, "job_scripts")
