FROM andreasga/cbctrecon:intel-opencl-runtime-cpu
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER root

# Useful for debugging docker: bash sudo vim \
RUN apt update \
    && apt install -y ninja-build clang cmake qt5-default \
        libfftw3-dev libtinfo-dev libeigen3-dev \
        libdcmtk-dev valgrind

RUN echo "export CC=clang" >> /etc/profile \
    && echo "export CXX=clang++" >> /etc/profile

USER user

RUN git config --global user.email "andreasga22@gmail.com" \
    && git config --global user.name "Andreas Gravgaard Andersen"

RUN echo "export CC=clang" >> ~/.bashrc \
    && echo "export CXX=clang++" >> ~/.bashrc

WORKDIR /home/user/
