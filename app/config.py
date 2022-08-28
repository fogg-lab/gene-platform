import os
from base64 import b64encode
import subprocess
from datetime import timedelta
from dotenv import load_dotenv

load_dotenv()

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
