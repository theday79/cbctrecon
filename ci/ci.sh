#/bin/bash
set -xe

if [[ "$CI_RUNNER_TAGS" = *"agravgaard-runner"* ]]; then
    MAKE_CONCURRENCY="-j13"
else
    MAKE_CONCURRENCY="-j2"
fi;

mkdir build && cd build

echo Test building: $BUILD_TESTING

# TinyRefl is for static reflection in the test, HDF5 is just to reduce complilation time:
if [ -d /home/user/ITK-build ]; then # Use system DCMTK, ITK and RTK
    cmake .. -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=ON -DUSE_TINYREFL=ON -DCBCTRECON_BUILD_TESTS=ON -DUSE_SYSTEM_DCMTK=ON -DUSE_SYSTEM_ITK=ON -DDCMTK_DIR=/home/user/DCMTK-build -DITK_DIR=/home/user/ITK-build -DRTK_DIR=/home/user/RTK-build
else
    cmake .. -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=OFF -DUSE_TINYREFL=ON -DCBCTRECON_BUILD_TESTS=ON -DITK_USE_SYSTEM_HDF5=ON
fi

cmake ..
make $MAKE_CONCURRENCY CbctRecon
cmake .. # Just to make sure tinyrefl data is generated
make $MAKE_CONCURRENCY CbctRecon_test
ctest -V
