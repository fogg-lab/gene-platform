FROM ubuntu:20.04
COPY install_packages.sh .
RUN ./install_packages.sh
