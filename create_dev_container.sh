#/bin/bash

# This script creates a development environment for the Singularity container.
# Prerequisite: Build image.sif, the Singularity image.

### Convert the production singularity container to a sandbox container
sudo singularity build --sandbox dev_container/ image.sif

### Install VSCode server in the sandbox container
sudo singularity exec --writable dev_container/ bash -c \
"curl -fsSL https://code-server.dev/install.sh | sh; \
ln -s /usr/local/bin/pip /usr/bin/pip3.8; \
pip install jupyterlab; \
pip install icecream;"

###############################################################################
# Using the development container
# 
# To start the container normally, run:
#   singularity run dev_container/
# 
# Or start the container in writable mode (e.g., to install other packages):
#   sudo singularity run --bind $PWD:/root --writable dev_container/
# 
# Example - Start the VSCode server in the container:
#   singularity run dev_container/
#   code-server --auth none
# Start the Jupyter server in the container:
#   singularity run dev_container/
#   jupyter-lab
###############################################################################
