#/bin/bash
set -xe

# Show available OpenCL driver & device:
clinfo

mkdir -p build && cd build

echo Test building: $BUILD_TESTING

export COMMON_FLAGS=".. -GNinja -DCMAKE_CXX_STANDARD=17 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/user/ -DBUILD_TESTING=OFF -DCBCTRECON_BUILD_TESTS=ON -DRTK_USE_OPENCL=ON"
# Eigen should be included in ITK if necessary:
export COMMON_SYSTEM_LIBS="-DUSE_SYSTEM_ZLIB=ON -DUSE_SYSTEM_dlib=ON -DUSE_SYSTEM_DCMTK=ON -DHUNTER_ENABLED=OFF -DUSE_SYSTEM_Plastimatch=OFF"

export FORWARD_FLAGS="-DCMAKE_CXX_FLAGS='$CXXFLAGS' -DCMAKE_EXE_LINKER_FLAGS='$LDFLAGS'"

if [ -d /home/user/ITK-build ]; then # Use system ITK and RTK
    export COMMON_NONSYSTEM_ITK="-DUSE_SYSTEM_ITK=ON -DITK_DIR=/home/user/ITK-build"
else
    # USE_SYSTEM_EIGEN = OFF with hunter disabled actually triggers the use of system eigen (will be corrected eventually)
    export COMMON_NONSYSTEM_ITK="-DUSE_SYSTEM_ITK=OFF -DITK_USE_SYSTEM_DCMTK=ON -DUSE_DCMTK=ON -DUSE_ITK_DCMTK=OFF -DUSE_SYSTEM_EIGEN=OFF"
fi

if [[ "$CUDA_AVAILABLE" = "YES" ]]; then
    export CUDA_FLAGS="-DUSE_CUDA=ON -DEXACT_GCC=/usr/bin/gcc-8"
    nvidia-smi
    # Ubuntu Bionic:
    export DLIBDIR="-Ddlib_DIR=/usr/lib/cmake/dlib/"
else
    export CUDA_FLAGS="-DUSE_CUDA=OFF"
    # Ubuntu Eoan:
    export DLIBDIR="-Ddlib_DIR=/usr/lib/x86_64-linux-gnu/cmake/dlib/"
fi

if [[ "$COVERAGE" = "YES" ]]; then
    export COVERAGE_FLAGS="-DCBCTRECON_COVERAGE=ON"
else
    export COVERAGE_FLAGS="-DCBCTRECON_COVERAGE=OFF"
fi

cmake $COMMON_FLAGS $FORWARD_FLAGS $COMMON_SYSTEM_LIBS $COMMON_NONSYSTEM_ITK $CUDA_FLAGS $COVERAGE_FLAGS $DLIBDIR
cmake --build . --target CbctData
cmake --build .

ctest -VV

if [[ "$COVERAGE" = "YES" ]]; then
    cmake --build . --target fastcov_html
#  lcov --directory Testing --base-directory ../Library/CbctReconLib/ --capture --output-file coverage.info
#  lcov --summary coverage.info
fi
