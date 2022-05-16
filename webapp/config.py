#Flask config

from datetime import timedelta

SESSION_PERMANENT = True
#SECRET_KEY = 'abc'
SESSION_TYPE = "filesystem"
SESSION_COOKIE_NAME = "my_session"
PERMANENT_SESSION_LIFETIME = timedelta(minutes=30)
debug = True
