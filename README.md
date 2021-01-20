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
 - [ITK](https://github.com/InsightSoftwareConsortium/ITK) (MSVC: Automatic download and configuration is broken, see "Build ITK for `USE_SYSTEM_ITK=ON`" below)
 - [RTK](https://github.com/SimonRit/RTK)
 - [Plastimatch](https://gitlab.com/plastimatch/plastimatch)
 - [dlib](https://github.com/davisking/dlib)

## goPMC Extension (The GPU monte carlo based proton dose calculation engine)
Before installation consider:
goPMC binary files and a Visual Studio 2013 Debug Compiler is needed for goPMC support
and the goPMC app must be compiled seperately.

## How to build
I recommend that you look at the short bash script, [ci/ci.sh](https://gitlab.com/agravgaard/cbctrecon/blob/master/ci/ci.sh), for an example of how to configure and build the project.

Possible config and build on linux, with CUDA:

```
git clone https://gitlab.com/agravgaard/cbctrecon.git
mkdir build
cd build

cmake -GNinja ../cbctrecon \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_INSTALL_PREFIX=/home/user/cbctrecon_install \
  -DHUNTER_ENABLED=ON \
  -DUSE_CUDA=ON \
  -DEXACT_GCC=/ust/bin/gcc-8 \
  -DUSE_SYSTEM_ITK=OFF \
  -DUSE_SYSTEM_DCMTK=OFF \
  -DUSE_SYSTEM_Plastimatch=OFF \
  -DUSE_SYSTEM_ZLIB=OFF \
  -DUSE_SYSTEM_dlib=OFF

ninja
ninja INSTALL
```

`EXACT_GCC` is needed for CUDA to work with newer versions of gcc or clang, it's an internal variable used in plastimatch and RTK, which may soon be deprecated.

`CMAKE_INSTALL_PREFIX` must be a path the user has write permissions in, otherwise temporary files can't be created at runtime.

(using ninja because it's much faster than Make and msbuild)

Possible config and build on Windows, with OpenCL 2.1:
In the "x64 Native Tools Command Prompt for VS 2019"

```
REM cd %userprofile% or whereever you want the project

git clone https://gitlab.com/agravgaard/cbctrecon.git
mkdir build
cd build

cmake -G"Visual Studio 16 2019" -A"x64" ../cbctrecon ^
  -DCMAKE_BUILD_TYPE=RelWithDebInfo ^
  -DCMAKE_PREFIX_PATH="C:/Qt/5.14.1/msvc2017_64/" ^
  -DCMAKE_INSTALL_PREFIX="C:/Users/user/cbctrecon_install/" ^
  -DHUNTER_ENABLED=ON ^
  -DRTK_USE_OPENCL=ON ^
  -DCBCTRECON_OPENCL_VERSION=210 ^
  -DITK_DIR="C:/Program Files (x86)/ITK/lib/cmake/ITK-5.1" ^
  -DUSE_SYSTEM_DCMTK=OFF ^
  -DUSE_SYSTEM_Plastimatch=OFF ^
  -DUSE_SYSTEM_ZLIB=OFF ^
  -DUSE_SYSTEM_dlib=OFF

cmake --build . --config RelWithDebInfo -j N
REM Where N is the number of CPU cores you want to assign to compiling

cmake --build . --config RelWithDebInfo --target INSTALL
```
Assuming git and cmake is in path

## Build ITK for `USE_SYSTEM_ITK=ON`

Certain option are necessary when you configure ITK manually:

```
// For RTK
Module_RTK=ON

// Optional for Cuda
Module_ITKCudaCommon=ON
RTK_USE_CUDA=ON

// For Plastimatch
Module_Review=ON
Module_ITKDeprecated=ON

// Optional but recommended for saving time and space:
BUILD_TESTING=OFF
ITK_BUILD_DEFAULT_MODULES=OFF
BUILD_EXAMPLES=OFF
ITK_USE_KWSTYLE=OFF

// Optional for more performance
ITK_USE_SYSTEM_FFTW=ON // Just download the binary, or use MKL or cuFFT
ITK_USE_FFTWD=ON
ITK_USE_FFTWF=ON
Module_TBB=ON // A bit difficult to configure, as you'll have to generate the TBBConfig.cmake yourself

```

## Build DCMTK for `USE_SYSTEM_DCMTK=ON`
```
DCMTK_OVERWRITE_WIN32i_COMPILER_FLAGS=OFF
```

## Build Plastimatch for `USE_SYSTEM_Plastimatch=ON`
These options are necessary to avoid linking errors (because the supertbuild would conflict with the ITK and DCMTK we need) and support C++ 17 (requires a newer version of dlib than the one included in plastimatch)
```
PLM_CONFIG_ENABLE_SUPERBUILD=OFF
ITK_DIR=/*Whereever you built or installed ITK*/
DCMTK_DIR=/*Whereever you built or installed DCMTK*/
dlib_DIR=/*wherever you installed dlib*/
```

## Avoiding DLL and linker hell
This is only a Windows problem. Make sure that all projects were compiled with the same linker option (shared / static) in all `CMAKE_CXX_FLAGS_*` and `CMAKE_C_FLAGS_*`, either `/MD` or `/MT`. These doesn't mix well. Prefer `/MT` and set `BUILD_SHARED_LIBS=OFF` in all projects, this seems to be the most stable configuration, although it's not forwarded by hunter (at least not to dlib).


