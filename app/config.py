import os
from datetime import timedelta
from dotenv import load_dotenv

load_dotenv()

SECRET_KEY = os.getenv("SESSION_SECRET_KEY", None)
if SECRET_KEY is None:
    os.environ["SESSION_SECRET_KEY"] = os.urandom(24)
    SECRET_KEY = os.getenv("SESSION_SECRET_KEY")

SESSION_PERMANENT = True
SESSION_TYPE = "filesystem"
SESSION_COOKIE_NAME = "my_session"
PERMANENT_SESSION_LIFETIME = timedelta(minutes=30)
debug = True
