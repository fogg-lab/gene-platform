import os.path
from sys import executable
from redis import Redis
from rq import Queue
from flask import current_app

from app.task_utils.task_runner import TaskRunner
from app.task_utils.execute_job_async import execute_job_async
from app.helper import StatusDict, get_valid_invalid_dsets


class PreprocessingRunner(TaskRunner):
    """Class for preparing and running a preprocessing task."""
    def __init__(self, task_id, task_dir):
        super().__init__(task_id, task_dir)
        self.task_type = "preprocessing"
        self._input_filenames = ["config.yml"]

    GEO_DATA_PATH = current_app.config["GEO_DATA_PATH"]
    GDC_DATA_PATH = current_app.config["GDC_DATA_PATH"]
    PROBEMAP_PATH = current_app.config["PROBEMAP_PATH"]
    SCRIPTS_PATH = current_app.config["SCRIPTS_PATH"]
    PREP_ALL_SCRIPT = os.path.join(SCRIPTS_PATH, "preprocess.py")
    PREP_GEO_SCRIPT = os.path.join(SCRIPTS_PATH, "prep_geo.r")
    PREP_GDC_SCRIPT = os.path.join(SCRIPTS_PATH, "prep_gdc.r")

    def execute_task(self) -> StatusDict:
        """Run a preprocessing task"""
        cfg = self.get_config()
        data_source = cfg.get("data_source")
        dsets = cfg.get("dsets")

        py_script = PreprocessingRunner.PREP_ALL_SCRIPT
        if "geo" in data_source.lower():
            r_script = PreprocessingRunner.PREP_GEO_SCRIPT
            data_dir = PreprocessingRunner.GEO_DATA_PATH
        elif "gdc" in data_source.lower():
            r_script = PreprocessingRunner.PREP_GDC_SCRIPT
            data_dir = PreprocessingRunner.GDC_DATA_PATH
        else:
            raise ValueError("Invalid data source")

        probemap_path = PreprocessingRunner.PROBEMAP_PATH
        input_dir = os.path.join(self._task_dir, "input")
        output_dir = os.path.join(self._task_dir, "output")
        log = os.path.join(self._task_dir, ".log")

        # Delete any previous preprocessing results
        for file in os.listdir(input_dir):
            if file.endswith("processed.tsv"):
                os.remove(os.path.join(input_dir, file))

        # Remake the config file with only valid datasets
        dsets, _ = get_valid_invalid_dsets(dsets, data_source, data_dir)
        cfg["dsets"] = dsets
        self.save_config(cfg)
        cfg_path = os.path.join(self._task_dir, "input/config.yml")

        cmd = (f"{executable} {py_script} -s 'Rscript {r_script}' -c {cfg_path} -d {data_dir} "
               f"-o {output_dir} -l {log}")

        if "geo" in data_source.lower():
            cmd += f" -p {probemap_path}"

        q = Queue(connection=Redis())
        q.enqueue(execute_job_async, args=(self._task_id, log, self.task_type, cmd),
                  job_timeout=600)

        return StatusDict(status="", errors=[])

    def validate_config(self, config) -> StatusDict:
        """Ensures config parameters are valid for the task"""
        source = config.get("data_source", "")
        dsets = config.get("dsets", "")

        if "gdc" in source.lower():
            metadata_dir = PreprocessingRunner.GDC_DATA_PATH
        elif "geo" in source.lower():
            metadata_dir = PreprocessingRunner.GEO_DATA_PATH
        else:
            metadata_dir = ""

        valid_dsets, invalid_dsets = get_valid_invalid_dsets(dsets, source, metadata_dir)

        return PreprocessingRunner._get_preprocessing_confirmation_msg(valid_dsets,
                                                                       invalid_dsets, source)

    def validate_task(self) -> StatusDict:
        """Validates all input files for the task"""
        return StatusDict(status="", errors=[])

    @staticmethod
    def _get_preprocessing_confirmation_msg(valid_dsets, invalid_dsets, source):
        """Return html confirmation message with datasets to load"""

        if "gdc" not in source.lower() and "geo" not in source.lower():
            err_msg = f"<p><b>Error:</b> Invalid source: {source}.</p>"
            return StatusDict(status=err_msg, errors=[])

        if not valid_dsets:
            err_msg = ("<p><b>Error:</b> No provided datasets were recognized from the source:"
                   f" {source}.</p>")
            return StatusDict(status=err_msg, errors=[])

        status_msg = ""

        if len(invalid_dsets) > 0:
            status_msg += "<p><b>Warning:</b> The following datasets were not recognized:"
            status_msg += "</b></p>"
            for dset in invalid_dsets:
                status_msg += f"<p>{dset}</p>"
            status_msg += "<br>"

        status_msg += f"<p><b>You are requesting the following datasets from {source}:"
        status_msg += "</b></p>"
        for dset in valid_dsets:
            status_msg += f"<p>{dset.upper()}</p>"
        status_msg += "<br>"
        status_msg += "<p>Proceed?</p>"

        return StatusDict(status=status_msg, errors=[])
