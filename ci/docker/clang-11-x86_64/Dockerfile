FROM andreasga/cbctrecon:intel-opencl-runtime-cpu
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER root

# Useful for debugging docker: bash sudo vim \
RUN apt update \
    && apt install -y \
        ninja-build g++ qt5-default \
        libfftw3-dev libtinfo-dev libeigen3-dev \
        libdcmtk-dev libgdcm2-dev valgrind \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN cd /usr/share/cmake/Modules/ \
    && cmake -D INSTALL_DIR=/usr/lib/cmake/TBB \
        -D SYSTEM_NAME=Linux \
        -D LIB_PATH=/usr/lib/x86_64-linux-gnu \
        -D INC_PATH=/usr/include \
        -P tbb_config_installer.cmake

# The flag -DLIBCXX_ENABLE_PARALLEL_ALGORITHMS=ON  creates problems in cmake config.
RUN git clone https://github.com/llvm/llvm-project.git --depth=1 --shallow-submodules \
    && cd llvm-project \
    && mkdir build \
    && cd build \
    && cmake -GNinja ../llvm \
        -D LLVM_ENABLE_PROJECTS="clang;lld;libcxx;libcxxabi;openmp;pstl" \
        -D LLVM_EXTERNAL_PSTL_SOURCE_DIR="../pstl" \
        -D CMAKE_BUILD_TYPE=Release \
        -D PSTL_PARALLEL_BACKEND="tbb" \
        -D LIBCXX_ENABLE_PARALLEL_ALGORITHMS=ON \
    && ninja \
    && ninja install \
    && cd ../.. \
    && rm -Rf ./llvm-project \
    && ldconfig \
    && /usr/local/bin/clang --version

RUN echo "export CC=/usr/local/bin/clang" >> /etc/profile \
    && echo "export CXX=/usr/local/bin/clang++" >> /etc/profile

USER user

RUN git config --global user.email "andreasga22@gmail.com" \
    && git config --global user.name "Andreas Gravgaard Andersen"

RUN echo "export CC=/usr/local/bin/clang" >> ~/.bashrc \
    && echo "export CXX=/usr/local/bin/clang++" >> ~/.bashrc

WORKDIR /home/user/
