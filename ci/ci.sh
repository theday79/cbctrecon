#/bin/bash
set -xe

# Show available OpenCL driver & device:
clinfo

mkdir -p build && cd build

echo Test building: $BUILD_TESTING

export COMMON_FLAGS=".. -GNinja -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=OFF -DCBCTRECON_BUILD_TESTS=ON"
# Eigen should be included in ITK if necessary:
export COMMON_SYSTEM_LIBS="-DUSE_SYSTEM_ZLIB=ON -DUSE_SYSTEM_DCMTK=ON -DHUNTER_ENABLED=OFF -DUSE_SYSTEM_Plastimatch=OFF"

if [ -d /home/user/ITK-build ]; then # Use system ITK and RTK
    export COMMON_NONSYSTEM_ITK="-DUSE_SYSTEM_ITK=ON -DITK_DIR=/home/user/ITK-build"
else
    # USE_SYSTEM_EIGEN = OFF with hunter disabled actually triggers the use of system eigen (will be corrected eventually)
    export COMMON_NONSYSTEM_ITK="-DUSE_SYSTEM_ITK=OFF -DITK_USE_SYSTEM_DCMTK=ON -DUSE_DCMTK=ON -DUSE_ITK_DCMTK=OFF -DUSE_SYSTEM_EIGEN=OFF"
fi

if [[ "$CUDA_AVAILABLE" = "YES" ]]; then
    export CUDA_FLAGS="-DUSE_CUDA=ON -DEXACT_GCC=/usr/bin/gcc-7"
    nvidia-smi
else
    export CUDA_FLAGS="-DUSE_CUDA=OFF"
fi

if [[ "$COVERAGE" = "YES" ]]; then
    export COVERAGE_FLAGS="-DCMAKE_BUILD_TYPE=Debug -DCBCTRECON_COVERAGE=ON"
else
    export COVERAGE_FLAGS="-DCMAKE_BUILD_TYPE=RelWithDebInfo -DCBCTRECON_COVERAGE=OFF"
fi

cmake $COMMON_FLAGS $COMMON_SYSTEM_LIBS $COMMON_NONSYSTEM_ITK $CUDA_FLAGS $COVERAGE_FLAGS
cmake --build . --target CbctData
cmake --build .

ctest -VV

if [[ "$COVERAGE" = "YES" ]]; then
  lcov --directory Testing --base-directory ../Library/CbctReconLib/ --capture --output-file coverage.info
  lcov --summary coverage.info
fi
