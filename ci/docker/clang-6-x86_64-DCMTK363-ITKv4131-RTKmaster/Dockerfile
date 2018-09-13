FROM ubuntu:18.04
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER root

# Useful for debugging docker: bash sudo vim \
RUN echo "deb http://us.archive.ubuntu.com/ubuntu/ bionic multiverse" >> /etc/apt/sources.list \
    && apt update \
    && apt install -y make clang git cmake qt5-default libfftw3-dev libsqlite3-dev ocl-icd-opencl-dev libtinfo-dev valgrind libhdf5-dev \ 
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
    && git config --global user.name "Andreas Gravgaard Andersen" \
    && export CC=clang && export CXX=clang++

# DCMTK
RUN cd /home/user/ \
    && git clone --depth 1 --branch "DCMTK-3.6.3" https://github.com/DCMTK/DCMTK.git DCMTK \
    && mkdir -p /home/user/DCMTK-build \
    && cd DCMTK-build \
    && cmake ../DCMTK \
      -DDCMTK_WITH_DOXYGEN=OFF -DDCMTK_ENABLE_CXX11=ON -DDCMTK_ENABLE_STL=ON \
      -DDCMTK_OVERWRITE_WIN32_COMPILER_FLAGS=OFF -DDCMTK_WITH_XML=OFF -DDCMTK_WITH_PNG=OFF \
      -DDCMTK_WITH_TIFF=OFF -DDCMTK_WITH_OPENSSL=OFF -DDCMTK_WITH_SNDFILE=OFF -DDCMTK_WITH_ICONV=OFF \
    && cmake --build . -- -j13

# ITK
RUN cd /home/user/ \
    && git clone --depth 1 --branch "v4.13.1" https://github.com/InsightSoftwareConsortium/ITK.git ITK \
    && mkdir -p /home/user/ITK-build \
    && cd ITK-build \
    && cmake ../ITK -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF \
      -DModule_ITKReview=ON -DITK_USE_GPU=ON -DITK_BUILD_DEFAULT_MODULES=OFF \
      -DITK_USE_SYSTEM_FFTW=ON -DITK_USE_FFTWF=ON -DITK_USE_FFTWD=ON \
      -DUSE_SYSTEM_DCMTK=ON -DDCMTK_DIR=/home/user/DCMTK-build -DITK_USE_SYSTEM_HDF5=ON \
      -DModule_ITKMathematicalMorphology=ON -DModule_ITKConvolution=ON -DModule_ITKIOCSV=ON -DModule_ITKImageNoise=ON \
    && cmake --build . -- -j13

# RTK
RUN cd /home/user/ \
    && git clone --depth 1 https://github.com/SimonRit/RTK.git RTK \
    && mkdir -p /home/user/RTK-build \
    && cd RTK-build \
    && cmake ../RTK -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF \
      -DRTK_USE_CUDA=OFF -DRTK_TIME_EACH_FILTER=ON -DRTK_BUILD_APPLICATIONS=OFF \
      -DITK_DIR=/home/user/ITK-build \
    && cmake --build . -- -j13


WORKDIR /home/user/
