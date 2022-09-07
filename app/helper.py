"""This module contains helper functions for the app"""

import os
import csv
import yaml


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


def dict_to_yaml(dictionary):
    """
    Converts a dictionary to bytes object ready to be written to a YML file.
    Args:
        dictionary (dict): Dictionary to convert.
    Returns:
        bytes: String ready to be written to a YML file.
    """

    return bytes(yaml.dump(dictionary), "utf-8")


def yaml_to_dict(yaml_filepath):
    """
    Converts a YAML file to a dictionary.
    Args:
        yaml_filepath (string): Path to a YAML file.
    Returns:
        dict: Dictionary representation of the YAML file.
    """

    with open(yaml_filepath, "r", encoding="UTF-8") as yaml_file:
        return yaml.safe_load(yaml_file)
