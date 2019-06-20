#/bin/bash
set -xe

# Show available OpenCL driver & device:
clinfo

mkdir -p build && cd build

echo Test building: $BUILD_TESTING

export COMMON_FLAGS=".. -GNinja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=OFF -DCBCTRECON_BUILD_TESTS=ON"
# Eigen should be included in ITK if necessary:
export COMMON_SYSTEM_LIBS="-DUSE_SYSTEM_ZLIB=ON -DHUNTER_ENABLED=OFF"

if [ -d /home/user/ITK-build ]; then # Use system DCMTK, ITK and RTK
    if [ -d /home/user/RTK-build ]; then
        cmake $COMMON_FLAGS \
            -DUSE_CUDA=OFF -DCBCTRECON_COVERAGE=OFF \
            -DUSE_SYSTEM_DCMTK=ON -DUSE_SYSTEM_ITK=ON -DUSE_SYSTEM_RTK=ON \
            -DDCMTK_DIR=/home/user/DCMTK-build -DITK_DIR=/home/user/ITK-build -DRTK_DIR=/home/user/RTK-build \
            $COMMON_SYSTEM_LIBS
    else
        cmake $COMMON_FLAGS \
            -DUSE_CUDA=OFF -DCBCTRECON_COVERAGE=OFF \
            -DCBCTRECON_BUILD_TESTS=ON -DUSE_SYSTEM_DCMTK=ON -DUSE_SYSTEM_ITK=ON \
            -DDCMTK_DIR=/home/user/DCMTK-build -DITK_DIR=/home/user/ITK-build \
            $COMMON_SYSTEM_LIBS
    fi
else
    export COMMON_NONSYSTEM_ITK="-DITK_USE_SYSTEM_DCMTK=ON -DUSE_DCMTK=ON -DUSE_ITK_DCMTK=OFF -DUSE_SYSTEM_EIGEN=ON"
    if [[ "$CUDA_AVAILABLE" = "YES" ]]; then
        if [[ "$COVERAGE" = "YES" ]]; then
            cmake $COMMON_FLAGS \
                -DUSE_CUDA=ON -DEXACT_GCC="/usr/bin/gcc-7" \
                -DCBCTRECON_COVERAGE=ON \
                $COMMON_SYSTEM_LIBS $COMMON_NONSYSTEM_ITK
        else
            cmake $COMMON_FLAGS
                -DUSE_CUDA=ON -DEXACT_GCC="/usr/bin/gcc-7" \
                -DCBCTRECON_COVERAGE=OFF \
                $COMMON_SYSTEM_LIBS $COMMON_NONSYSTEM_ITK
        fi
    else
        cmake $COMMON_FLAGS \
            -DUSE_CUDA=OFF -DCBCTRECON_COVERAGE=OFF \
            $COMMON_SYSTEM_LIBS $COMMON_NONSYSTEM_ITK
    fi
fi

cmake ..
cmake --build . --target CbctData
cmake --build .

ctest -VV

if [[ "$COVERAGE" = "YES" ]]; then
  lcov --directory Testing --base-directory ../Library/CbctReconLib/ --capture --output-file coverage.info
  lcov --summary coverage.info
fi
