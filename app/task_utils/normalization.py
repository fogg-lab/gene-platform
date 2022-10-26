"""Class for preparing and running a normalization task."""
from app.task_utils.task_runner import TaskRunner


class NormalizationRunner(TaskRunner):
    def __init__(self, task_id, task_dir):
        super().__init__(task_id, task_dir)
        self.task_type = "normalization"
        self._input_filenames = ["counts.tsv", "coldata.tsv", "config.yml"]


    def update_task(self):
        """Task has a new input file - perform input validation."""
        pass


    def start_task(self):
        """Run a normalization task"""
        pass
