FROM gliderlabs/alpine:edge
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

# build-base depends on: gcc g++ make libc-dev fortify-headers
RUN apk update && \
    apk upgrade && \
    apk --update add --no-cache build-base libstdc++ bash git fftw-dev qt5-qtbase-dev cmake opencl-headers opencl-icd-loader opencl-icd-loader-dev && \
    git config --global user.email "andreasga22@gmail.com" && \
    git config --global user.name "Andreas Gravgaard Andersen" && \
    mkdir -p /home/user/cbctrecon/build && \
    git clone https://gitlab.com/agravgaard/cbctrecon.git /home/user/cbctrecon && \
    cd /home/user/cbctrecon && \
    git checkout ModernizeCMake && \
    git pull

USER user
ENV HOME=/home/user
WORKDIR /home/user/cbctrecon/build
