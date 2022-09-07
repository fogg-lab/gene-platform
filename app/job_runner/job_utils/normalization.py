"""
Functions for preparing and running normalization.
Used by the job runner module.
"""
import os
import subprocess
import csv
import time
import yaml
from flask import current_app

INPUT_FNAMES = ["counts.tsv", "coldata.tsv", "config.yml"]


def update_job(directory):
    """Job has a new input file - perform input validation."""
    pass


def start_job(directory):
    """Run a normalization job"""
    pass
