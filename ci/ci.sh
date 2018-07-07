#/bin/bash
set -xe

if [[ "$CI_RUNNER_TAGS" = *"agravgaard-runner"* ]]; then
    MAKE_CONCURRENCY="-j13"
else
    MAKE_CONCURRENCY="-j3"
fi;

mkdir build && cd build

echo Test building: $BUILD_TESTING

if [ "$BUILD_TESTING" == "YES" ]; then
    cmake .. -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=ON
    cmake .. 
    make $MAKE_CONCURRENCY
else
    cmake .. -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=OFF
    cmake .. 
    make $MAKE_CONCURRENCY
fi
