#/bin/bash
set -xe

if [[ "$CI_RUNNER_TAGS" = *"agravgaard-runner"* ]]; then
    MAKE_CONCURRENCY="-j13"
else
    MAKE_CONCURRENCY="-j2"
fi;

# Show available OpenCL driver & device:
clinfo

mkdir build && cd build

echo Test building: $BUILD_TESTING

# TinyRefl is for static reflection in the test, HDF5 is just to reduce complilation time:
if [ -d /home/user/ITK-build ]; then # Use system DCMTK, ITK and RTK
    if [ -d /home/user/RTK-build ]; then
        cmake .. -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=ON -DCBCTRECON_BUILD_TESTS=ON -DUSE_SYSTEM_DCMTK=ON -DUSE_SYSTEM_ITK=ON -DUSE_SYSTEM_RTK=ON -DDCMTK_DIR=/home/user/DCMTK-build -DITK_DIR=/home/user/ITK-build -DRTK_DIR=/home/user/RTK-build
    else
        cmake .. -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=ON -DCBCTRECON_BUILD_TESTS=ON -DUSE_SYSTEM_DCMTK=ON -DUSE_SYSTEM_ITK=ON -DDCMTK_DIR=/home/user/DCMTK-build -DITK_DIR=/home/user/ITK-build
    fi
else
    if [[ "$CUDA_AVAILABLE" = "YES" ]]; then
        if [[ "$COVERAGE" = "YES" ]]; then
            cmake .. -DUSE_CUDA=ON -DEXACT_GCC="/usr/bin/gcc-7" -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=OFF -DCBCTRECON_BUILD_TESTS=ON -DITK_USE_SYSTEM_HDF5=ON -DCMAKE_BUILD_TYPE=Debug -DCBCTRECON_COVERAGE=ON
        else
            cmake .. -DUSE_CUDA=ON -DEXACT_GCC="/usr/bin/gcc-7" -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=OFF -DCBCTRECON_BUILD_TESTS=ON -DITK_USE_SYSTEM_HDF5=ON
        fi
    else
        cmake .. -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=OFF -DCBCTRECON_BUILD_TESTS=ON -DITK_USE_SYSTEM_HDF5=ON
    fi
fi

cmake ..
make $MAKE_CONCURRENCY CbctRecon
cmake ..
make $MAKE_CONCURRENCY CbctData
make $MAKE_CONCURRENCY CbctRecon_test
ctest -V
if [[ "$COVERAGE" = "YES" ]]; then
    make $MAKE_CONCURRENCY CbctReconLib_coverage
fi
