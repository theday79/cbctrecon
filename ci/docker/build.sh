#!/bin/bash

REPO="andreasga"
NAME="cbctrecon"

CC="gcc-7-x86_64"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

DCMIRTK="DCMTK363-ITKv4131-RTKmaster"
docker build --tag ${REPO}/${NAME}:${CC}-${DCMIRTK} docker/${CC}-${DCMIRTK} 
docker push ${REPO}/${NAME}:${CC}-${DCMIRTK} 

CC="clang-6-x86_64"
docker build --tag ${REPO}/${NAME}:${CC} docker/${CC}
docker push ${REPO}/${NAME}:${CC}

docker build --tag ${REPO}/${NAME}:${CC}-${DCMIRTK} docker/${CC}-${DCMIRTK} 
docker push ${REPO}/${NAME}:${CC}-${DCMIRTK} 

# docker run -it --rm ${REPO}/${NAME}:gcc7-x86_64 bash
