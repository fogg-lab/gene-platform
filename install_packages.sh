#!/bin/bash

set -euo pipefail

export DEBIAN_FRONTEND=noninteractive
apt -y update
apt -y upgrade

export TZ=Etc/UTC
export LC_ALL=C.UTF-8

apt install -y dirmngr gnupg apt-transport-https ca-certificates \
    software-properties-common libxml2-dev libssl-dev libcurl4-openssl-dev \
    cmake git libfontconfig1-dev xclip libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev curl

apt-key adv --keyserver keyserver.ubuntu.com --recv-keys \
    E298A3A825C0D65DFD57CBB651716619E084DAB9
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
curl -fsSL https://packages.redis.io/gpg | gpg --dearmor -o \
    /usr/share/keyrings/redis-archive-keyring.gpg
echo "deb [signed-by=/usr/share/keyrings/redis-archive-keyring.gpg] \
    https://packages.redis.io/deb $(lsb_release -cs) main" | tee \
    /etc/apt/sources.list.d/redis.list

add-apt-repository ppa:deadsnakes/ppa

apt -y update

apt install -y r-base redis python3.11 python3.11-distutils python3.11-dev
curl -sS https://bootstrap.pypa.io/get-pip.py | python3.11

python3.11 -m pip install --upgrade pip
python3.11 -m pip install --upgrade setuptools
python3.11 -m pip install Flask Flask-Session Flask-Login Flask-SQLAlchemy pyyaml \
    numpy pandas gunicorn oauthlib requests python-dotenv[cli] rq

R -e "install.packages(c('yaml', 'ggpubr', 'tidyverse', 'corrplot', 'corrr', \
    'BiocManager'), repos='https://cloud.r-project.org/')"
R -e "BiocManager::install(c('BiocParallel', 'remotes', 'BioinformaticsFMRP/TCGAbiolinks', \
    'HDF5Array', 'affycoretools', 'affyPLM', 'sva', 'edgeR', 'GEOquery', 'biomaRt'))"
R -e "install.packages('devtools', repos='https://cloud.r-project.org/')"
R -e "devtools::install_github('zhangyuqing/sva-devel')"

apt-get clean

rm -rf /var/lib/apt/lists/*

unset DEBIAN_FRONTEND
