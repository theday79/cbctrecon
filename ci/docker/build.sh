#!/bin/bash

REPO="andreasga"
NAME="cbctrecon"

CC="conan-python"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

CC="intel-opencl-runtime-cpu"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

CC="gcc-9-x86_64"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

DCMIRTK="libDCMTK-ITK-RTK"
echo "docker build --tag ${REPO}/${NAME}:${CC}-${DCMIRTK} docker/${CC}-${DCMIRTK}"
docker build --tag ${REPO}/${NAME}:${CC}-${DCMIRTK} docker/${CC}-${DCMIRTK}
docker push ${REPO}/${NAME}:${CC}-${DCMIRTK}

CC="gcc-9-CUDA-x86_64"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

CC="clang-11-x86_64"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

echo "docker build --tag ${REPO}/${NAME}:${CC}-${DCMIRTK} docker/${CC}-${DCMIRTK}"
docker build --tag ${REPO}/${NAME}:${CC}-${DCMIRTK} docker/${CC}-${DCMIRTK}
docker push ${REPO}/${NAME}:${CC}-${DCMIRTK}

CC="clang-11-ROCm-x86_64"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

CC="icc-21-oneAPI-x86_64"
echo "docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

# docker run -it --rm ${REPO}/${NAME}:gcc7-x86_64 bash
