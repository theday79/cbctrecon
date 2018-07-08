#!/bin/bash

REPO="andreasga"
NAME="cbctrecon"

docker build --tag ${REPO}/${NAME}:gcc7-x86_64 docker/gcc-7-x86_64
docker push ${REPO}/${NAME}:gcc7-x86_64

# docker run -it --rm ${REPO}/${NAME}:gcc7-x86_64 bash
