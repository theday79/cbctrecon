#!/bin/bash

REPO="andreasga"
NAME="cbctrecon"

CC="conan-python"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC} || { echo "failed" ; exit 1; }

CC="intel-opencl-runtime-cpu"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC} || { echo "failed" ; exit 1; }

CC="gcc-10-x86_64"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC} || { echo "failed" ; exit 1; }

DCMIRTK="libDCMTK-ITK-RTK"
echo "docker build --tag ${REPO}/${NAME}:${CC}-${DCMIRTK} docker/${CC}-${DCMIRTK}"
docker build --tag ${REPO}/${NAME}:${CC}-${DCMIRTK} docker/${CC}-${DCMIRTK}
docker push ${REPO}/${NAME}:${CC}-${DCMIRTK} || { echo "failed" ; exit 1; }

## commented out until RTK fixes compute flags for CUDA 11
CC="gcc-10-CUDA-x86_64"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC} || { echo "failed" ; exit 1; }

CC="clang-11-x86_64"
if [ "$1" != "" ]; then
    echo "LLVM will be recompiled"
    echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
    docker build --tag ${REPO}/${NAME}:${CC} docker/${CC} || { echo "failed" ; exit 1; }
    docker push ${REPO}/${NAME}:${CC} || { echo "failed" ; exit 1; }
else
    echo "LLVM will NOT be recompiled, use 'sudo ./docker/build.sh build_llvm' if you want to recompile"
fi

echo "docker build --tag ${REPO}/${NAME}:${CC}-${DCMIRTK} docker/${CC}-${DCMIRTK}"
docker build --tag ${REPO}/${NAME}:${CC}-${DCMIRTK} docker/${CC}-${DCMIRTK}
docker push ${REPO}/${NAME}:${CC}-${DCMIRTK} || { echo "failed" ; exit 1; }

CC="clang-11-ROCm-x86_64"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC} || { echo "failed" ; exit 1; }

CC="icc-21-oneAPI-x86_64"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC} || { echo "failed" ; exit 1; }

