#/bin/bash
set -xe

if [[ "$CI_RUNNER_TAGS" = *"agravgaard-runner"* ]]; then
    MAKE_CONCURRENCY="-j13"
else
    MAKE_CONCURRENCY="-j2"
fi;

mkdir build && cd build

echo Test building: $BUILD_TESTING

if [ "$BUILD_TESTING" == "YES" ]; then
    cmake .. -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=ON -DUSE_TINYREFL=ON -DCBCTRECON_BUILD_TESTS=ON
    cmake .. 
    make $MAKE_CONCURRENCY CbctRecon CbctRecon_test
else
    cmake .. -DUSE_CUDA=OFF -DCMAKE_INSTALL_PREFIX="/home/user/" -DBUILD_TESTING=OFF -DUSE_TINYREFL=ON -DCBCTRECON_BUILD_TESTS=ON
    cmake .. 
    make $MAKE_CONCURRENCY CbctRecon CbctRecon_test
fi
