FROM ubuntu:18.04
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER root

# Useful for debugging docker: bash sudo vim \
RUN echo "deb http://us.archive.ubuntu.com/ubuntu/ bionic multiverse" >> /etc/apt/sources.list \
    && apt update \
    && apt install -y make gcc g++ git cmake qt5-default libfftw3-dev libsqlite3-dev ocl-icd-opencl-dev libtinfo-dev valgrind libhdf5-dev clinfo \ 
    && groupadd 1001 -g 1001 \
    && groupadd 1000 -g 1000 \
    && groupadd 2000 -g 2000 \
    && groupadd 999 -g 999 \
    && useradd -ms /bin/bash user -g 1001 -G 1000,2000,999 \
    && printf "user:user" | chpasswd \
    && adduser user sudo \
    && printf "user ALL= NOPASSWD: ALL\\n" >> /etc/sudoers

RUN mkdir neo && cd neo \
    && wget https://github.com/intel/compute-runtime/releases/download/18.47.11882/intel-gmmlib_18.4.0.348_amd64.deb \
    && wget https://github.com/intel/compute-runtime/releases/download/18.47.11882/intel-igc-core_18.44.1060_amd64.deb \
    && wget https://github.com/intel/compute-runtime/releases/download/18.47.11882/intel-igc-opencl_18.44.1060_amd64.deb \
    && wget https://github.com/intel/compute-runtime/releases/download/18.47.11882/intel-opencl_18.47.11882_amd64.deb \
    && dpkg -i *.deb

USER user

RUN git config --global user.email "andreasga22@gmail.com" \
    && git config --global user.name "Andreas Gravgaard Andersen"

WORKDIR /home/user/
