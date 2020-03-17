FROM andreasga/cbctrecon:intel-opencl-runtime-cpu
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER root

# Useful for debugging docker: bash sudo vim \
RUN apt update \
    && apt install -y ninja-build qt5-default \
        libfftw3-dev libtinfo-dev libeigen3-dev \
        libdcmtk-dev valgrind

RUN git clone https://github.com/llvm/llvm-project.git --depth=1 --shallow-submodules \
    && cd llvm-project \
    && mkdir build \
    && cd build \
    && cmake -GNinja ../llvm \
        -DLLVM_ENABLE_PROJECTS="clang;libcxx;libcxxabi;openmp;pstl" \
        -DCMAKE_BUILD_TYPE=Release \
    && ninja \
    && ninja install \
    && /usr/local/bin/clang --version

RUN echo "export CC=/usr/local/bin/clang" >> /etc/profile \
    && echo "export CXX=/usr/local/bin/clang++" >> /etc/profile

USER user

RUN git config --global user.email "andreasga22@gmail.com" \
    && git config --global user.name "Andreas Gravgaard Andersen"

RUN echo "export CC=/usr/local/bin/clang" >> ~/.bashrc \
    && echo "export CXX=/usr/local/bin/clang++" >> ~/.bashrc

WORKDIR /home/user/
