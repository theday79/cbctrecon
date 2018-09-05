#!/bin/bash

REPO="andreasga"
NAME="cbctrecon"

docker build --tag ${REPO}/${NAME}:gcc7-x86_64 docker/gcc-7-x86_64
docker push ${REPO}/${NAME}:gcc7-x86_64

docker build --tag ${REPO}/${NAME}:gcc7-x86_64-DCMTK3.6.3-ITKv4.13.1-RTKmaster docker/gcc-7-x86_64-DCMTK3.6.3-ITKv4.13.1-RTKmaster 
docker push ${REPO}/${NAME}:gcc7-x86_64-DCMTK3.6.3-ITKv4.13.1-RTKmaster 

# docker run -it --rm ${REPO}/${NAME}:gcc7-x86_64 bash
