# GENE Platform

## Gene Expression Explorer (GENE) Platform: A web app for exploring gene expression data and running differential expression analyses.

### Main site (currently you must sign in with your OSU account): https://geneplatform.oregonstate.edu

## Installation
The web server runs on Linux, including the Windows Subsystem for Linux (WSL).

## Installation

Prerequisite: Install [Docker](https://docs.docker.com/get-docker/) or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) container platform

Instructions to set up and start the web app:
1. Clone this repo
```bash
git clone https://github.com/fogg-lab/gene-platform.git
```

2. Build the container
```bash
cd gene-platform
chmod +x install_packages.sh
# One of the following commands (depends on whether you are using Singularity or Docker):
sudo singularity build gene-platform-image.sif singularity.def
# or, for docker:
docker build -t gene-platform-image .
```

3. Run the container
```bash
# One of the following:
singularity run gene-platform-image.sif
# or, for docker:
# Replace 1234 with some port number (the same number should be used in the next step)
docker run -it --volume $PWD:/gene-platform -w /gene-platform -p 1234:1234 gene-platform-image
```

4. Run the web app (replace 1234 with some port number)
```bash
gunicorn --workers=3 --threads=2 --bind 0.0.0.0:1234 wsgi:app
# ...or run the server persistently in the background like this:
nohup gunicorn --workers=3 --threads=2 --bind 0.0.0.0:1234 wsgi:app &
```
