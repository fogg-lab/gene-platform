class TaskNotFound(Exception):
    """Raised when a user-associated task is not found in the database"""

class InvalidTaskType(Exception):
    """Raised when a task type is not valid"""

class InvalidTaskInputFile(Exception):
    """Raised when an input filename is not recognized for the task type"""

class InvalidTaskOutput(Exception):
    """Raised when an output file is not recognized for the task type"""
