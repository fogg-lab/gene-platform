"""This module contains helper functions for tasks and flask blueprints."""

import os
import csv
from pathlib import Path
from typing import TypedDict, List
import json


def get_tsv_rows(filepath: str) -> List[List[str]]:
    """
    Returns the rows of the user input file as a 2d array
    If the file does not exist, returns None
    Args:
        filepath (string): Path to a .tsv file.
    Returns:
        List[List[str]]: Contents of the file as a 2d array.
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
        msg_list (List[str]): List of status messages.
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
    errors: List[str]


def get_series_probesets(accessions, metadata_dir):
    """Returns a dict of accessions with their respective probesets"""

    with open(Path(metadata_dir) / "series_platforms.json") as series_platforms_json:
        series_platforms = json.load(series_platforms_json)

    with open(Path(metadata_dir) / "platform_probesets.json") as platform_probesets_json:
        platform_probesets = json.load(platform_probesets_json)

    supported_probesets = set()
    for probemap_fname in [Path(f) for f in os.listdir(Path(metadata_dir) / "probe_maps")]:
        fname_base = probemap_fname.stem
        probeset = fname_base.split("_", 1)[1]
        supported_probesets.add(probeset)

    series_probesets = {}

    for accession in accessions:
        accession = accession.lower()
        if accession in series_platforms.keys():
            platform = series_platforms[accession]
        else:
            continue
        if (platform in platform_probesets.keys() and
                platform_probesets[platform] in supported_probesets):
            series_probesets[accession] = platform_probesets[platform]
        for probeset in platform_probesets.values():
            if platform in probeset and probeset in supported_probesets:
                series_probesets[accession] = probeset

    return series_probesets


def get_valid_geo_accessions(accessions, metadata_dir):
    """
    Returns a list of valid GEO accessions from given list of accessions
    Args:
        accessions (List[str]): List of GEO accessions
    Returns:
        list: List of valid GEO accessions
    """

    valid_accessions = list(get_series_probesets(accessions, metadata_dir).keys())

    return valid_accessions


def get_valid_gdc_projects(project_names, metadata_dir):
    """Returns a list of valid GDC projects from project_names"""

    with open(Path(metadata_dir) / "gdc_projects.json") as gdc_projects_json:
        gdc_projects = json.load(gdc_projects_json)

    valid_projects = []
    for project_name in project_names:
        if "-" not in project_name:
            continue
        proj_name_lower = project_name.lower()
        proj_base, proj_ext = proj_name_lower.split("-", 1)
        if proj_base in gdc_projects and proj_ext in gdc_projects[proj_base]:
            valid_projects.append(project_name)

    return valid_projects


def get_valid_invalid_dsets(dsets, source, metadata_dir):
    """Return valid and invalid datasets from user input"""
    if " " in dsets and not "," in dsets:
        dsets = dsets.split(" ")
    elif "," in dsets:
        dsets = dsets.split(",")
    else:
        dsets = [dsets.lower()]
    for i, dset in enumerate(dsets):
        dsets[i] = dset.replace(" ", "").lower()

    if "gdc" in source.lower():
        valid_dsets = get_valid_gdc_projects(dsets, metadata_dir)
    elif "geo" in source.lower():
        valid_dsets = get_valid_geo_accessions(dsets, metadata_dir)

    print(type(dsets))
    print(dsets)
    print(len(dsets))
    print(type(valid_dsets))
    print(valid_dsets)
    print(len(valid_dsets))

    invalid_dsets = list(set(dsets).difference(set(valid_dsets)))

    return valid_dsets, invalid_dsets
