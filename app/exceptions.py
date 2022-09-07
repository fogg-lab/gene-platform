class JobNotFoundError(Exception):
    """Raised when a user-associated job is not found in the database"""

class InvalidJobType(Exception):
    """Raised when a job type is not valid"""

class InvalidInputFile(Exception):
    """Raised when an input filename is not recognized for the job type"""
