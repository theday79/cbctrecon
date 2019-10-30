FROM andreasga/cbctrecon:gcc-7-x86_64
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER root

RUN apt-get update \
    && apt-get install -y software-properties-common \
    && add-apt-repository ppa:ubuntu-toolchain-r/test \
    && apt-get update \
    && apt-get install -y gcc-9 g++-9 lcov

USER user

# ITK
RUN export CC=/usr/bin/gcc-9 \
    && export CXX=/usr/bin/g++-9 \
    && cd /home/user/ \
    && git clone --depth 1 --branch "v5.0.1" https://github.com/InsightSoftwareConsortium/ITK.git ITK \
    && mkdir -p /home/user/ITK-build \
    && cd ITK-build \
    && cmake -GNinja ../ITK -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF \
      -DModule_ITKReview=ON -DITK_USE_GPU=ON -DITK_BUILD_DEFAULT_MODULES=OFF \
      -DITK_USE_SYSTEM_FFTW=ON -DITK_USE_FFTWF=ON -DITK_USE_FFTWD=ON \
      -DModule_ITKDCMTK=ON -DITK_USE_SYSTEM_DCMTK=ON -DModule_ITKDeprecated=ON \
      -DITK_USE_SYSTEM_ZLIB=ON -DModule_RTK=ON \
      -DRTK_USE_CUDA=OFF -DRTK_TIME_EACH_FILTER=ON -DRTK_BUILD_APPLICATIONS=OFF \
      -DREMOTE_GIT_TAG_RTK="master" \
    && cmake --build .

WORKDIR /home/user/
