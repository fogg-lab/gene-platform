# syntax=docker/dockerfile:1
FROM r-base:latest
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends build-essential libpq-dev python3.8 python3-pip python3-setuptools python3-dev
RUN pip3 install --upgrade pip

ENV PYTHONPATH "${PYTHONPATH}:/app"
WORKDIR /app

ADD requirements.txt .
ADD requirements.r .

# installing python libraries
RUN pip3 install -r requirements.txt

# installing r libraries
RUN Rscript requirements.r
