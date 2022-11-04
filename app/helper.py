"""This module contains helper functions for tasks and flask blueprints."""

import os
import csv
from typing import TypedDict


def get_tsv_rows(filepath):
    """
    Returns the rows of the user input file as a 2d array
    If the file does not exist, returns None
    Args:
        filepath (string): Path to a .tsv file.
    Returns:
        list[list[str]]: Contents of the file as a 2d array.
    """

    rows = None

    if os.path.isfile(filepath):
        with open(filepath, "r", encoding="UTF-8") as tsv_file:
            rows = list(csv.reader(tsv_file, delimiter="\t"))

    return rows

def format_msg_list_html(msg_list, tag="p"):
    """
    Formats a list of status messages into HTML for display in the browser.
    Args:
        msg_list (list[str]): List of status messages.
        tag (str): HTML tag to use for each message (li or p). Defaults to p.
    Returns:
        str: HTML-formatted status messages.
    """

    if tag not in ("li", "p"):
        raise ValueError("tag must be li or p")

    if msg_list is None or len(msg_list) == 0:
        return ""

    opening_tag = f"<{tag}>"
    closing_tag = f"</{tag}>"

    msg_html = "".join([f"{opening_tag}{msg}{closing_tag}" for msg in msg_list])

    return msg_html

def get_analysis_confirmation_msg(config_params):
    """Get analysis formula from the config parameters"""

    contrast_level = config_params["contrast_level"]
    reference_level = config_params["reference_level"]

    analysis_formula =  "<p><b>You are performing the following analysis:</b></p>\n"
    analysis_formula += f"<p><i>condition ~ (intercept) + {contrast_level}</i></p>\n\n"
    analysis_formula += f"<p>where the reference group is <i>{reference_level}</i>.\n</p>"

    return analysis_formula

class StatusDict(TypedDict):
    """Status dictionary for task operations that need to return a status message."""
    status: str
    errors: list[str]
