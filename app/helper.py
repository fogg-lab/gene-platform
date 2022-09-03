"""
This module contains miscellaneous utility functions for the app
"""

import os
import csv


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
