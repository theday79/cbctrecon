FROM ubuntu:18.04
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER root

# Useful for debugging docker: bash sudo vim \
RUN echo "deb http://us.archive.ubuntu.com/ubuntu/ bionic multiverse" >> /etc/apt/sources.list \
    && apt update \
    && apt install -y make gcc g++ git cmake qt5-default libfftw3-dev libsqlite3-dev ocl-icd-opencl-dev \ 
    && groupadd 1001 -g 1001 \
    && groupadd 1000 -g 1000 \
    && groupadd 2000 -g 2000 \
    && groupadd 999 -g 999 \
    && useradd -ms /bin/bash user -g 1001 -G 1000,2000,999 \
    && printf "user:user" | chpasswd \
    && adduser user sudo \
    && printf "user ALL= NOPASSWD: ALL\\n" >> /etc/sudoers


USER user

RUN git config --global user.email "andreasga22@gmail.com" \
    && git config --global user.name "Andreas Gravgaard Andersen" 
#    && mkdir -p /home/user/cbctrecon \
#    && git clone https://gitlab.com/agravgaard/cbctrecon.git /home/user/cbctrecon \
#    && cd /home/user/cbctrecon \
#    && git checkout ModernizeCMake \
#    && git pull \
#    && mkdir build \
#    && cd build \
#    && cmake .. -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user" \
#    && cmake .. 


WORKDIR /home/user/cbctrecon
