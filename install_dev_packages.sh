#!/bin/bash

# This script creates a Singularity container with packages for development.
# Prerequisite: Built image.sif, the Singularity image for production.
# Usage: ./create_dev_container.sh

# 1: Convert the production singularity container to a sandbox container
sudo rm -rf tmp_container/
sudo singularity build --sandbox tmp_container/ image.sif

# 2: Install development packages in the sandbox container
sudo singularity exec --writable tmp_container/ bash -c \
"apt -y purge python3.8; \
python3.11 -m pip install executing; \
python3.11 -m pip install icecream; \
python3.11 -m pip install jupyterlab; 
curl -fsSL https://code-server.dev/install.sh | sh;" \

# 3: Convert the sandbox container to an sif file
rm -rf dev_image.sif
sudo singularity build dev_image.sif tmp_container/

# 4: Delete the temporary sandbox container
sudo rm -rf tmp_container/
