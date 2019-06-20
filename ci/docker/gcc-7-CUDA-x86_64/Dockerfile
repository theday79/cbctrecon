# Latest = 10.0 devel ubuntu 18.04, devel is FROM runtime which is FROM base (30/11/2018)
FROM nvidia/cuda:latest
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER root

## OpenCL driver (With modification from nvidia/opencl)
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ocl-icd-libopencl1 \
        ocl-icd-opencl-dev \
        clinfo

RUN mkdir -p /etc/OpenCL/vendors \
    && echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd

# nvidia-container-runtime
ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES compute,utility
## end OpenCL driver

# Useful for debugging docker: bash sudo vim \
RUN echo "deb http://us.archive.ubuntu.com/ubuntu/ bionic multiverse" >> /etc/apt/sources.list \
    && apt update \
    && apt install -y --no-install-recommends \
        ninja-build gcc g++ git wget qt5-default libqt5opengl5-dev \
        libfftw3-dev libsqlite3-dev libtinfo-dev libhdf5-dev \
        libeigen3-dev libdcmtk-dev zlib1g-dev valgrind python lcov \
    && rm -rf /var/lib/apt/lists/*

# Note: we specify an exact cmake version to work around a cmake 3.10 CUDA 10 issue.
# Reference: https://github.com/clab/dynet/issues/1457
RUN mkdir /opt/cmake && cd /opt/cmake \
    && wget -nv https://cmake.org/files/v3.14/cmake-3.14.5-Linux-x86_64.sh \
    && sh cmake-3.14.5-Linux-x86_64.sh --prefix=/opt/cmake --skip-license \
    && ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake \
    && ln -s /opt/cmake/bin/ctest /usr/local/bin/ctest \
    && rm cmake-3.14.5-Linux-x86_64.sh \
    && cmake --version

RUN groupadd 1001 -g 1001 \
    && groupadd 1000 -g 1000 \
    && groupadd 2000 -g 2000 \
    && groupadd 999 -g 999 \
    && useradd -ms /bin/bash user -g 1001 -G 1000,2000,999 \
    && printf "user:user" | chpasswd \
    && adduser user sudo \
    && printf "user ALL= NOPASSWD: ALL\\n" >> /etc/sudoers


USER user
# Plastimatch won't download unless a git user is given:
RUN git config --global user.email "andreasga22@gmail.com" \
    && git config --global user.name "Andreas Gravgaard Andersen"

WORKDIR /home/user/
