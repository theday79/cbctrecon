FROM andreasga/cbctrecon:clang-11-x86_64
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER root

# Register the ROCM package repository, and install rocm-dev package
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends curl libnuma-dev gnupg \
    && curl -sL http://repo.radeon.com/rocm/apt/debian/rocm.gpg.key | apt-key add - \
    && printf "deb [arch=amd64] http://repo.radeon.com/rocm/apt/debian/ xenial main" | tee /etc/apt/sources.list.d/rocm.list \
    && apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        libelf1 rocm-dev build-essential \
        qt5-default lcov openssl libgdbm-dev \
        libgdbm-compat-dev liblzma-dev libffi-dev uuid-dev \
        software-properties-common apt-utils \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN cd /opt/rocm \
    && ln -s /usr/local/ llvm

# Build ROCm device library
RUN export PATH=/opt/rocm/llvm/bin/:$PATH \
    && git clone -b amd-stg-open https://github.com/RadeonOpenCompute/ROCm-Device-Libs.git \
    && cd ROCm-Device-Libs \
    && mkdir -p build && cd build \
    && CC=clang CXX=clang++ cmake -G Ninja .. \
        -D LLVM_DIR=/opt/rocm/llvm/ \
        -D CMAKE_BUILD_TYPE=Release \
        -D LLVM_ENABLE_ASSERTIONS=1 \
    && ninja \
    && ninja install

# Build HIP (not hip-hcc)
RUN apt update \
    && apt install -y --no-install-recommends libtinfo5 kmod \
    && rm -rf /var/lib/apt/lists/* \
    && git clone -b master https://github.com/ROCm-Developer-Tools/HIP.git \
    && cd HIP \
    && mkdir build && cd build \
    && cmake -G Ninja .. \
        -D CMAKE_INSTALL_PREFIX=/opt/rocm/hip \
        -D HIP_COMPILER=clang \
        -D CMAKE_BUILD_TYPE=Release \
    && ninja \
    && ninja install \
    && cd ../.. \
    && rm -Rf ./HIP

ENV HCC_AMDGPU_TARGET=gfx900

# Verify installation:
RUN /opt/rocm/hip/bin/hipconfig --full

# hipSYCL
RUN git clone https://github.com/illuhad/hipSYCL --recurse-submodules \
  && cd hipSYCL \
  && mkdir build && cd build \
  && cmake ../ -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DLLVM_DIR=/usr/local/lib/cmake/llvm/ \
    -DDISABLE_LLVM_VERSION_CHECK=ON \
    -DWITH_ROCM_BACKEND=ON \
    -DHIPCC_COMPILER=/opt/rocm/bin/hipcc \
    -DCLANG_INCLUDE_PATH=/usr/local/lib/clang/11.0.0/ \
    -DWITH_CPU_BACKEND=ON \
  && ninja \
  && ninja install \
  && cd ../.. \
  && rm -Rf hipSYCL

ENV HIPSYCL_GPU_ARCH=gfx900

# SyclParallelSTL
RUN git clone https://github.com/KhronosGroup/SyclParallelSTL.git \
  && cd SyclParallelSTL \
  && sed -i '/endif/a target_include_directories(${ADD_SYCL_TARGET} PUBLIC /usr/local/include/hipSYCL/contrib/)' /usr/local/lib/cmake/hipsycl-config.cmake \
  && sed -i 's/include(FindTriSYCL)/find_package(hipsycl REQUIRED)/g' ./CMakeLists.txt \
  && sed -i 's/triSYCL/hipSYCL/g' ./CMakeLists.txt \
  && sed -i 's/sycl_to_target(/sycl_to_target(TARGET /g' ./examples/CMakeLists.txt \
  && sed -i 's/sycl_to_target(/sycl_to_target(TARGET /g' ./tests/pstl-tests/CMakeLists.txt \
  && mkdir build && cd build \
  && cmake ../ -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_CXX_FLAGS="-fopenmp" \
    -DUSE_COMPUTECPP=OFF \
    -DHIPSYCL_PLATFORM=rocm \
    -DHIPSYCL_GPU_ARCH=gfx900 \
    -Dhipsycl_DIR=/usr/local/lib/cmake/ \
  && ninja \
  && ninja install \
  && cd ../.. \
  && rm -Rf SyclParallelSTL


USER user

# Plastimatch won't download unless a git user is given:
RUN git config --global user.email "andreasga22@gmail.com" \
    && git config --global user.name "Andreas Gravgaard Andersen"

# Build ITK and RTK by themselves without CUDA support, although we may be able to hipify it, I doubt it's worth the effort
# "master" because 5.1rc01 are missing some set/get functions required for cuda part of RTK
# Also, we'll have to use my patched fork of RTK to forward the correct std=c++ flag to CUDA
RUN cd /home/user/ \
    && git clone --depth 1 --branch "master" https://github.com/InsightSoftwareConsortium/ITK.git ITK \
    && sed -i 's/SimonRit/agravgaard/g' ./ITK/Modules/Remote/RTK.remote.cmake \
    && sed -i 's/GIT_TAG .*/GIT_TAG master/g' ./ITK/Modules/Remote/RTK.remote.cmake \
    && sed -i 's/<<</<< </g' ./ITK/Modules/Core/Common/include/itkArray.h \
    && sed -i 's/<<</<< </g' ./ITK/Modules/Core/Common/src/itkArrayOutputSpecialization.cxx \
    && sed -i 's/defined(H5_ZLIB_HEADER)/true/g' ./ITK/Modules/ThirdParty/HDF5/src/itkhdf5/src/H5Zdeflate.c \
    && sed -i 's/H5_ZLIB_HEADER/\"itk_zlib.h\"/g' ./ITK/Modules/ThirdParty/HDF5/src/itkhdf5/src/H5Zdeflate.c \
    && mkdir -p /home/user/ITK-build \
    && cd ITK-build \
    && export CC=/opt/rocm/hip/bin/hipcc  \
    && export CXX=/opt/rocm/hip/bin/hipcc \
    && cmake -G Ninja ../ITK \
      -D CMAKE_CXX_STANDARD=17 \
      -D CMAKE_BUILD_TYPE=Release \
      -D BUILD_TESTING=OFF -D BUILD_EXAMPLES=OFF \
      -D Module_ITKReview=ON \
      -D ITK_BUILD_DEFAULT_MODULES=OFF \
      -D ITK_USE_SYSTEM_FFTW=ON -D ITK_USE_FFTWF=ON -D ITK_USE_FFTWD=ON \
      -D Module_ITKDCMTK=ON \
      -D ITK_USE_SYSTEM_DCMTK=ON \
      -D ITK_USE_SYSTEM_GDCM=ON \
      -D Module_ITKDeprecated=ON \
      -D ITK_USE_SYSTEM_ZLIB=ON \
      -D Module_RTK=ON \
      -D Module_ITKTBB=ON \
      -D RTK_USE_CUDA=OFF \
      -D RTK_BUILD_APPLICATIONS=OFF \
      -D ITK_MINIMUM_COMPLIANCE_LEVEL=3 \
    && ninja \
    && sudo ninja install \
    && cd ../ \
    && rm -Rf ITK*

WORKDIR /home/user/
