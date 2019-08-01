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
    && echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd \
    && ln -s /usr/local/cuda/lib64/stubs/libcuda.so /usr/local/cuda/lib64/stubs/libcuda.so.1 \
    && LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64/stubs/ \
    && export LD_LIBRARY_PATH \
    && ldconfig

# nvidia-container-runtime
ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES compute,utility
## end OpenCL driver

# Useful for debugging docker: bash sudo vim \
RUN echo "deb http://us.archive.ubuntu.com/ubuntu/ bionic multiverse" >> /etc/apt/sources.list \
    && apt update \
    && apt install -y ninja-build gcc g++ git qt5-default \
        libfftw3-dev libtinfo-dev libeigen3-dev \
        libdcmtk-dev valgrind zlib1g-dev wget clinfo lcov \
        sudo libsqlite3-dev sqlite3 bzip2 libbz2-dev libssl-dev openssl libgdbm-dev \
        libgdbm-compat-dev liblzma-dev libreadline-dev libncursesw5-dev libffi-dev uuid-dev \
    && rm -rf /var/lib/apt/lists/*
# From libsqlite3-dev and forward is just explicit pip dependencies

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

# Conani:
ENV PYENV_ROOT=/opt/pyenv \
    PATH=/opt/pyenv/shims:${PATH}

RUN wget --no-check-certificate --quiet -O /tmp/pyenv-installer https://github.com/pyenv/pyenv-installer/raw/master/bin/pyenv-installer \
    && chmod +x /tmp/pyenv-installer \
    && /tmp/pyenv-installer \
    && rm /tmp/pyenv-installer \
    && update-alternatives --install /usr/bin/pyenv pyenv /opt/pyenv/bin/pyenv 100 \
    && PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.7.1 \
    && pyenv global 3.7.1 \
    && pip install -q --upgrade --no-cache-dir pip \
    && pip install -q --no-cache-dir conan conan-package-tools \
    && chown -R user:1001 /opt/pyenv \
    # remove all __pycache__ directories created by pyenv
    && find /opt/pyenv -iname __pycache__ -print0 | xargs -0 rm -rf \
    && update-alternatives --install /usr/bin/python python /opt/pyenv/shims/python 100 \
    && update-alternatives --install /usr/bin/python3 python3 /opt/pyenv/shims/python3 100 \
    && update-alternatives --install /usr/bin/pip pip /opt/pyenv/shims/pip 100 \
    && update-alternatives --install /usr/bin/pip3 pip3 /opt/pyenv/shims/pip3 100


USER user

RUN mkdir -p /home/user/.conan \
    && printf 'eval "$(pyenv init -)"\n' >> ~/.bashrc \
    && printf 'eval "$(pyenv virtualenv-init -)"\n' >> ~/.bashrc

# Plastimatch won't download unless a git user is given:
RUN git config --global user.email "andreasga22@gmail.com" \
    && git config --global user.name "Andreas Gravgaard Andersen"

# Build ITK and RTK by themselves with CUDA support
RUN cd /home/user/ \
    && git clone --depth 1 --branch "v5.0.0" https://github.com/InsightSoftwareConsortium/ITK.git ITK \
    && mkdir -p /home/user/ITK-build \
    && cd ITK-build \
    && cmake -GNinja ../ITK -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF \
      -DModule_ITKReview=ON -DITK_USE_GPU=ON -DITK_BUILD_DEFAULT_MODULES=OFF \
      -DITK_USE_SYSTEM_FFTW=ON -DITK_USE_FFTWF=ON -DITK_USE_FFTWD=ON \
      -DModule_ITKDCMTK=ON -DITK_USE_SYSTEM_DCMTK=ON -DModule_ITKDeprecated=ON \
      -DITK_USE_SYSTEM_ZLIB=ON -DModule_RTK=ON -DModule_ITKCudaCommon=ON \
      -DRTK_USE_CUDA=ON -DRTK_TIME_EACH_FILTER=ON -DRTK_BUILD_APPLICATIONS=OFF \
      -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda -DCMAKE_LIBRARY_PATH=/usr/local/cuda/lib64/stubs/ \
    && cmake --build .

WORKDIR /home/user/
