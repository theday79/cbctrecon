FROM andreasga/cbctrecon:clang-6-x86_64
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER user

# DCMTK
RUN cd /home/user/ \
    && git clone --depth 1 --branch "DCMTK-3.6.3" https://github.com/DCMTK/DCMTK.git DCMTK \
    && mkdir -p /home/user/DCMTK-build \
    && cd DCMTK-build \
    && cmake -GNinja ../DCMTK \
      -DDCMTK_WITH_DOXYGEN=OFF -DDCMTK_ENABLE_CXX11=ON -DDCMTK_ENABLE_STL=ON \
      -DDCMTK_OVERWRITE_WIN32_COMPILER_FLAGS=OFF -DDCMTK_WITH_XML=OFF -DDCMTK_WITH_PNG=OFF \
      -DDCMTK_WITH_TIFF=OFF -DDCMTK_WITH_OPENSSL=OFF -DDCMTK_WITH_SNDFILE=OFF -DDCMTK_WITH_ICONV=OFF \
    && cmake --build . -- -j13

# ITK
RUN cd /home/user/ \
    && git clone --depth 1 --branch "v4.13.1" https://github.com/InsightSoftwareConsortium/ITK.git ITK \
    && echo "itk_fetch_module(RTK \"RTK external module for ITK\" GIT_REPOSITORY \"https://github.com/SimonRit/RTK.git\" GIT_TAG master)" > ITK/Modules/Remote/ExternalRTK.remote.cmake \
    && mkdir -p /home/user/ITK-build \
    && cd ITK-build \
    && cmake -GNinja ../ITK -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF \
      -DModule_ITKReview=ON -DITK_USE_GPU=ON -DITK_BUILD_DEFAULT_MODULES=OFF \
      -DITK_USE_SYSTEM_FFTW=ON -DITK_USE_FFTWF=ON -DITK_USE_FFTWD=ON \
      -DUSE_SYSTEM_DCMTK=ON -DDCMTK_DIR=/home/user/DCMTK-build -DITK_USE_SYSTEM_HDF5=ON \
      -DModule_ITKMathematicalMorphology=ON -DModule_ITKConvolution=ON -DModule_ITKIOCSV=ON -DModule_ITKImageNoise=ON \
      -DModule_RTK=ON -DRTK_USE_CUDA=OFF -DRTK_TIME_EACH_FILTER=ON -DRTK_BUILD_APPLICATIONS=OFF \
    && cmake --build . -- -j13


WORKDIR /home/user/
