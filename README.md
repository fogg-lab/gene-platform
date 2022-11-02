# GENE Platform

## Gene Expression Explorer (GENE) Platform: A web app for exploring gene expression data and running differential expression analyses.

### Main site (currently you must sign in with your OSU account): https://geneplatform.oregonstate.edu

## Installation

Prerequisite: Install [Docker](https://docs.docker.com/get-docker/) or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) container platform

Instructions to set up and start the web app:  
1. Clone this repo  
  `git clone git@github.com:wigginno/gene-platform.git`  

2. Build the container  
  `cd gene-platform`  
  `chmod +x install_packages.sh`  
  `\# One of the following commands (depends on whether you are using Docker or Singularity):`  
  `docker build -t gene-platform-image .`  
  `sudo singularity build gene-platform-image.sif singularity.def`  
  *Note*: It can take an hour or longer to build the image.  

3. Run the container  
  \# One of the following:  
  `singularity run gene-platform-image.sif`  
  `docker run gene-platform-image`

4. Run the web app (replace xxxx with a port number)  
  `gunicorn --workers=3 --threads=2 --bind 0.0.0.0:xxxx wsgi:app`  
  ...or run the server persistently in the background like this:  
  `nohup gunicorn --workers=3 --threads=2 --bind 0.0.0.0:xxxx wsgi:app &`  
  Note: Sometimes it may take 2 tries to start the server.  
