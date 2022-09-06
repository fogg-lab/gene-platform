import os
import json
import sys
import requests
from flask import Blueprint, request, redirect, current_app
from oauthlib.oauth2 import WebApplicationClient
from flask_login import login_user, logout_user, login_required, current_user
from dotenv import load_dotenv
from app.models.user import User

auth_bp = Blueprint('auth_bp', __name__)

load_dotenv(os.path.join(current_app.root_path, ".env"))
client_id = os.environ.get("GOOGLE_CLIENT_ID", None)
client_secret = os.environ.get("GOOGLE_CLIENT_SECRET", None)
redirect_uri = os.environ.get("REDIRECT_URI", None)

client = WebApplicationClient(client_id)


@auth_bp.before_app_request
def before_request():
    """Check if user is logged in on each page load"""
    is_login_endpoint = request.endpoint in ["auth_bp.login", "auth_bp.callback"]
    if not current_user.is_authenticated and not is_login_endpoint:
        print("redirecting...")
        return redirect("/login")


@auth_bp.route("/login")
def login():
    """Login page"""
    # Find out what URL to hit for Google login
    google_provider_cfg = get_google_provider_cfg()
    authorization_endpoint = google_provider_cfg["authorization_endpoint"]

    # Use library to construct the request for Google login and provide
    # scopes that let you retrieve user's profile from Google
    load_dotenv("../.env")
    print(f"REDIRECT_URI: {redirect_uri}", file=sys.stderr)
    print(f"authorization_endpoint: {authorization_endpoint}")
    request_uri = client.prepare_request_uri(
        authorization_endpoint,
        redirect_uri = redirect_uri,
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
    token_url, headers, body = client.prepare_token_request(
        token_endpoint,
        authorization_response=request.url.replace("http:","https:"),
        redirect_url=request.base_url.replace("http:","https:"),
        code=code
    )

    token_response = requests.post(
        token_url,
        headers=headers,
        data=body,
        auth=(client_id, client_secret),
    )

    # Parse the tokens
    client.parse_request_body_response(json.dumps(token_response.json()))

    # Get user info
    userinfo_endpoint = google_provider_cfg["userinfo_endpoint"]
    uri, headers, body = client.add_token(userinfo_endpoint)
    userinfo_response = requests.get(uri, headers=headers, data=body)

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
    login_user(user)

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
    return requests.get(discover_url).json()
