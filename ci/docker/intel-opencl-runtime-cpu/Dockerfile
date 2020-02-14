FROM andreasga/cbctrecon:conan-python
MAINTAINER Andreas Gravgaard Andersen <andreasga22@gmail.com>

USER root

RUN apt update \
    && echo "Install dependencies for Intel OpenCL runtime" \
    && apt install -y numactl libtbb2 zlib1g libncurses5 wget tar gzip bsdtar \
        ocl-icd-opencl-dev clinfo libtbb-dev libdlib-dev

RUN echo "Install Intel runtime for Intel Core and Xeon processors (Inspired by the AUR package intel-opencl-runtime)" \
    && wget "http://registrationcenter-download.intel.com/akdlm/irc_nas/13793/l_opencl_p_18.1.0.013.tgz" \
    && echo "208806279b0b9219ca6a17c64cbe0e4a3876a8b5d3f172bf296d85c0f1c74126  l_opencl_p_18.1.0.013.tgz" > test.sha256 \
    && sha256sum -c test.sha256 \
    && tar -xf l_opencl_p_18.1.0.013.tgz \
    && cd l_opencl_p_18.1.0.013 \
    && install -Dm644 license.txt /usr/share/licenses/intel-opencl-runtime/license \
    && rm rpm/intel-openclrt-pset-*.rpm \
    && for i in rpm/*.rpm; do bsdtar -xf "$i"; done \
    && mkdir -p /etc/OpenCL/vendors \
    && echo "/opt/intel/opencl-runtime/linux/compiler/lib/intel64_lin/libintelocl.so" > /etc/OpenCL/vendors/intel.icd \
    && mkdir -p /opt/intel/opencl-runtime \
    && cp -r opt/intel/opencl_*/* /opt/intel/opencl-runtime \
    && cd ..

