# CbctRecon

CBCT Reconstruction toolkit for Elekta and Varian type cone beam projections

[![pipeline status](https://gitlab.com/agravgaard/cbctrecon/badges/master/pipeline.svg)](https://gitlab.com/agravgaard/cbctrecon/commits/master)
[![coverage report](https://gitlab.com/agravgaard/cbctrecon/badges/master/coverage.svg)](https://gitlab.com/agravgaard/cbctrecon/commits/master)

## Proton dose calculation on scatter-corrected CBCT image: Feasibility study for adaptive proton therapy
http://dx.doi.org/10.1118/1.4923179

## Supported OS/compilers
 - Windows 10
   - MSVC 2019
 - Linux Manjaro rolling/Ubuntu Eoan and MacOS High Sierra
   - GCC (8) 9+
   - Clang 9.0.1+

Some compilers and OSs are not tested as often.
I expect the code to compile on most systems with a C++17 compatible compiler.

(clang 9.0.0 has has an ICE when compiling Abseil charconv, which we need, it's fixed in 9.0.1)

In order to compile the software, you must have installed the following prerequisites.

## Prerequisites
 - Git
 - Cmake >3.14
 - Qt 5.X
   - On Windows Set `CMAKE_PREFIX_PATH="[Qt install dir]/5.[X.Y]/msvc20[ZZ]_64/"`
 - OpenCL 
   - Use `-DRTK_USE_OPENCL=ON` to get OpenCL versions of the CUDA algorithms

## Optional, but recommended, dependencies
 - FFTW - Add directory to PATH
 - CUDA > 8.0
   - C++ compiler compatible with your version of CUDA. (See [this table](https://gist.github.com/ax3l/9489132))

## Below is downloaded and compiled automatically, if `USE_SYSTEM_XXXX=OFF`
 - [DCMTK](https://github.com/DCMTK/DCMTK)
 - [ITK](https://github.com/InsightSoftwareConsortium/ITK)
 - [RTK](https://github.com/SimonRit/RTK)
 - [Plastimatch](https://gitlab.com/plastimatch/plastimatch)
 - [dlib](https://github.com/davisking/dlib)

## goPMC Extension (The GPU monte carlo based proton dose calculation engine)
Before installation consider:
goPMC binary files and a Visual Studio 2013 Debug Compiler is needed for goPMC support
and the goPMC app must be compiled seperately.

## How to build
I recommend that you look at the short bash script, [ci/ci.sh](https://gitlab.com/agravgaard/cbctrecon/blob/master/ci/ci.sh), for an example of how to configure and build the project.

Possible config and build (on linux, with CUDA):

```
git clone https://gitlab.com/agravgaard/cbctrecon.git
mkdir build
cd build

cmake -GNinja ../cbctrecon \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DHUNTER_ENABLED=ON \
  -DUSE_CUDA=ON \
  -DEXACT_GCC=/ust/bin/gcc-8 \
  -DUSE_SYSTEM_ITK=OFF \
  -DUSE_SYSTEM_DCMTK=OFF \
  -DUSE_SYSTEM_Plastimatch=OFF \
  -DUSE_SYSTEM_ZLIB=OFF \
  -DUSE_SYSTEM_dlib=OFF

ninja
```

`EXACT_GCC` is needed for CUDA to work with newer versions of gcc or clang, it's an internal variable used in plastimatch and RTK, which may soon be deprecated.

(using ninja because it's much faster than Make and msbuild)
