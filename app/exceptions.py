class JobNotFound(Exception):
    """Raised when a user-associated job is not found in the database"""

class InvalidJobType(Exception):
    """Raised when a job type is not valid"""

class InvalidJobInputFile(Exception):
    """Raised when an input filename is not recognized for the job type"""

class InvalidJobOutput(Exception):
    """Raised when an output file is not recognized for the job type"""
