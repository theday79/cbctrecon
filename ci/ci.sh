#/bin/bash
set -xe

# Show available OpenCL driver & device:
clinfo

mkdir -p build && cd build

echo Test building: $BUILD_TESTING

# TinyRefl is for static reflection in the test, HDF5 is just to reduce complilation time:
if [ -d /home/user/ITK-build ]; then # Use system DCMTK, ITK and RTK
    if [ -d /home/user/RTK-build ]; then
        cmake .. -GNinja -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=ON -DCBCTRECON_BUILD_TESTS=ON -DUSE_SYSTEM_DCMTK=ON -DUSE_SYSTEM_ITK=ON -DUSE_SYSTEM_RTK=ON -DDCMTK_DIR=/home/user/DCMTK-build -DITK_DIR=/home/user/ITK-build -DRTK_DIR=/home/user/RTK-build
    else
        cmake .. -GNinja -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=ON -DCBCTRECON_BUILD_TESTS=ON -DUSE_SYSTEM_DCMTK=ON -DUSE_SYSTEM_ITK=ON -DDCMTK_DIR=/home/user/DCMTK-build -DITK_DIR=/home/user/ITK-build
    fi
else
    if [[ "$CUDA_AVAILABLE" = "YES" ]]; then
        if [[ "$COVERAGE" = "YES" ]]; then
            cmake .. -GNinja -DUSE_CUDA=ON -DEXACT_GCC="/usr/bin/gcc-7" -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=OFF -DCBCTRECON_BUILD_TESTS=ON -DITK_USE_SYSTEM_HDF5=ON -DCMAKE_BUILD_TYPE=Debug -DCBCTRECON_COVERAGE=ON
        else
            cmake .. -GNinja -DUSE_CUDA=ON -DEXACT_GCC="/usr/bin/gcc-7" -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=OFF -DCBCTRECON_BUILD_TESTS=ON -DITK_USE_SYSTEM_HDF5=ON
        fi
    else
        cmake .. -GNinja -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=OFF -DCBCTRECON_BUILD_TESTS=ON -DITK_USE_SYSTEM_HDF5=ON
    fi
fi

cmake ..
cmake --build . --target CbctData
cmake --build .

ctest -VV

if [[ "$COVERAGE" = "YES" ]]; then
    cmake --build . --target CbctReconLib_coverage
fi
