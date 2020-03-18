FROM andreasga/cbctrecon:intel-opencl-runtime-cpu
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER root

# Useful for debugging docker: bash sudo vim \
RUN apt update \
    && apt install -y g++-9 gcc-9 ninja-build qt5-default \
        libfftw3-dev libtinfo-dev \
        libdcmtk-dev libeigen3-dev valgrind

USER user

RUN git config --global user.email "andreasga22@gmail.com" \
    && git config --global user.name "Andreas Gravgaard Andersen"

WORKDIR /home/user/
