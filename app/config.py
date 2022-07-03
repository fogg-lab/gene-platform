#Flask config

from datetime import timedelta
import os

SECRET_key = ""
if os.path.isfile("session_secret_key.txt"):
    with open("session_secret_key.txt", "r") as f:
        SECRET_KEY = f.read()
else:
    with open("session_secret_key.txt", "w") as f:
        SECRET_KEY = os.urandom(16)
        f.write(str(SECRET_KEY))

SESSION_PERMANENT = True
SESSION_TYPE = "filesystem"
SESSION_COOKIE_NAME = "my_session"
PERMANENT_SESSION_LIFETIME = timedelta(minutes=30)
debug = True
