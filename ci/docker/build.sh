#!/bin/bash

REPO="andreasga"
NAME="cbctrecon"

CC="conan-python"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

CC="intel-opencl-runtime-cpu"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

CC="gcc-9-x86_64"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

DCMIRTK="libDCMTK-ITKv5-RTK"
docker build --tag ${REPO}/${NAME}:${CC}-${DCMIRTK} docker/${CC}-${DCMIRTK}
docker push ${REPO}/${NAME}:${CC}-${DCMIRTK}

CC="gcc-8-CUDA-x86_64"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

CC="clang-9-x86_64"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

docker build --tag ${REPO}/${NAME}:${CC}-${DCMIRTK} docker/${CC}-${DCMIRTK}
docker push ${REPO}/${NAME}:${CC}-${DCMIRTK}

# docker run -it --rm ${REPO}/${NAME}:gcc7-x86_64 bash
