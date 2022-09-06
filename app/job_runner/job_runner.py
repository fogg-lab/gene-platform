from app.job_runner.job_utils import (analysis, batch_correction, correlation,
                                      normalization, preprocessing,
                                      expression_data_validation)

def run_job(job_dir, job_type):
    """Start job"""
    if job_type == "analysis":
        analysis.run_job(job_dir)
    elif job_type == "batch_correction":
        batch_correction.run_job(job_dir)
    elif job_type == "correlation":
        correlation.run_job(job_dir)
    elif job_type == "normalization":
        normalization.run_job(job_dir)
    elif job_type == "preprocessing":
        preprocessing.run_job(job_dir)
    else:
        raise Exception("Invalid job type")

def update_job(job_dir, job_type):
    """
    Post an update to a job to trigger new input validation
    Return status message
    """
    if job_type == "analysis":
        return analysis.update_job(job_dir)
    elif job_type == "batch_correction":
        return batch_correction.update_job(job_dir)
    elif job_type == "correlation":
        return correlation.update_job(job_dir)
    elif job_type == "normalization":
        return normalization.update_job(job_dir)
    elif job_type == "preprocessing":
        return preprocessing.update_job(job_dir)
    else:
        raise Exception("Invalid job type")
