from app.task_utils.task_runner import TaskRunner


class NormalizationRunner(TaskRunner):
    """Class for preparing and running a normalization task."""
    def __init__(self, task_id, task_dir):
        super().__init__(task_id, task_dir)
        self.task_type = "normalization"
        self._input_filenames = ["counts.tsv", "coldata.tsv", "config.yml"]

    def execute_task(self):
        """Run a normalization task"""
        return dict(status="", warnings=[], errors=[])

    def validate_config(self, config) -> dict:
        """Ensures config parameters are valid for the task"""
        return dict(status="", warnings=[], errors=[])

    def validate_task(self) -> dict:
        """Validates all input files for the task"""
        return dict(status="", warnings=[], errors=[])
