import os
import json
import requests
from flask import Blueprint, request, redirect, current_app
from oauthlib.oauth2 import WebApplicationClient
from flask_login import login_user, logout_user, login_required, current_user
from dotenv import load_dotenv
from app.models.user import User

auth_bp = Blueprint('auth_bp', __name__)

auth_bp.client_id=None
auth_bp.client_secret=None
auth_bp.redirect_uri=None
auth_bp.client=None
auth_bp.initialized=False


def init_env():
    """Load environment variables"""
    load_dotenv(os.path.join(current_app.root_path, ".env"))
    auth_bp.client_id = os.environ.get("GOOGLE_CLIENT_ID", None)
    auth_bp.client_secret = os.environ.get("GOOGLE_CLIENT_SECRET", None)
    auth_bp.redirect_uri = os.environ.get("REDIRECT_URI", None)
    auth_bp.client = WebApplicationClient(auth_bp.client_id)
    auth_bp.initialized = True


@auth_bp.before_app_request
def before_request():
    """Check if user is logged in on each page load"""
    if not auth_bp.initialized:
        init_env()
    enable_google_auth = os.getenv("ENABLE_GOOGLE_AUTH")
    if enable_google_auth and enable_google_auth.lower() not in ["false", "0"]:
        is_login_endpoint = request.endpoint in ["auth_bp.login", "auth_bp.callback"]
        if not current_user.is_authenticated and not is_login_endpoint:
            return redirect("/login")
        if not is_login_endpoint:
            # Log the user in as a guest
            if not User.get("guest1"):
                User.create("guest1", "John Doe", "email@domain.com")
            user = User.get("guest1")
            login_user(user)


@auth_bp.route("/login")
def login():
    """Login page"""
    # Find out what URL to hit for Google login
    google_provider_cfg = get_google_provider_cfg()
    authorization_endpoint = google_provider_cfg["authorization_endpoint"]

    # Use library to construct the request for Google login and provide
    # scopes that let you retrieve user's profile from Google
    load_dotenv("../.env")
    request_uri = auth_bp.client.prepare_request_uri(
        authorization_endpoint,
        redirect_uri = auth_bp.redirect_uri,
        scope=["openid", "email", "profile"],
    )
    return redirect(request_uri)


@auth_bp.route("/login/callback")
def callback():
    """Callback function for Google OAuth"""

    # Get authorization code Google sent back to you
    code = request.args.get("code")

    # Get the token endpoint URL
    google_provider_cfg = get_google_provider_cfg()
    token_endpoint = google_provider_cfg["token_endpoint"]

    # Prepare and send a request to get tokens
    token_url, headers, body = auth_bp.client.prepare_token_request(
        token_endpoint,
        authorization_response=request.url.replace("http:","https:"),
        redirect_url=request.base_url.replace("http:","https:"),
        code=code
    )

    token_response = requests.post(
        token_url,
        headers=headers,
        data=body,
        auth=(auth_bp.client_id, auth_bp.client_secret),
        timeout=10,
    )

    # Parse the tokens
    auth_bp.client.parse_request_body_response(json.dumps(token_response.json()))

    # Get user info
    userinfo_endpoint = google_provider_cfg["userinfo_endpoint"]
    uri, headers, body = auth_bp.client.add_token(userinfo_endpoint)
    userinfo_response = requests.get(uri, headers=headers, data=body, timeout=10)

    if userinfo_response.json().get("email_verified"):
        unique_id = userinfo_response.json()["sub"]
        users_email = userinfo_response.json()["email"]
        users_name = userinfo_response.json()["given_name"]
    else:
        return "User email not available or not verified by Google.", 400

    # Create a user in your db with the information provided
    # by Google
    user = User(
        id_=unique_id, name=users_name, email=users_email
    )

    # Doesn't exist? Add it to the database.
    if not User.get(unique_id):
        User.create(unique_id, users_name, users_email)

    # Begin user session by logging the user in
    login_user(user, remember=True)

    # Send user back to homepage
    return redirect("/")


@auth_bp.route("/logout")
@login_required
def logout():
    """Logout current user"""
    logout_user()
    return redirect("/")


def get_google_provider_cfg():
    """Get Google provider configuration"""
    discover_url = os.getenv("GOOGLE_DISCOVERY_URL")
    return requests.get(discover_url, timeout=10).json()
