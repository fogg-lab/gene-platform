"""Class for preparing and running a normalization job."""
from app.job_utils.job_runner import JobRunner


class NormalizationRunner(JobRunner):
    def __init__(self, job_id, job_dir):
        super().__init__(job_id, job_dir)
        self.job_type = "normalization"
        self._input_filenames = ["counts.tsv", "coldata.tsv", "config.yml"]


    def update_job(self):
        """Job has a new input file - perform input validation."""
        pass


    def start_job(self):
        """Run a normalization job"""
        pass
