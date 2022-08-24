# GENE Platform

## Gene Expression Explorer (GENE) Platform: A web app for exploring gene expression data and running differential expression analyses.

## Installation
The web server runs on Linux, including the Windows Subsystem for Linux (WSL).

Prerequisite - Singularity container platform (click [here](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) for installation instructions).

Instructions to set up and start the web app:  
1. Clone this repo  
  `git clone git@github.com:wigginno/gene-platform.git`  

2. Build the container  
  `cd gene-platform`  
  `singularity build image.sif singularity_builder.def`  
  Note: It may take an hour or longer to build the image.  

3. Run the container  
  `singularity run image.sif`  

4. Run the web app (replace xxxx with a port number)  
  `gunicorn --workers=3 --threads=2 --bind 0.0.0.0:xxxx wsgi:app`  
  ...or run the server persistently in the background like this:  
  `nohup gunicorn --workers=3 --threads=2 --bind 0.0.0.0:xxxx wsgi:app &`  
  Note: Sometimes it may take 2 tries to start the server.  
