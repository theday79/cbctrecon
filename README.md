# CbctRecon
CBCT Reconstruction toolkit for Elekta and Varian type cone beam projections

[![pipeline status](https://gitlab.com/agravgaard/cbctrecon/badges/master/pipeline.svg)](https://gitlab.com/agravgaard/cbctrecon/commits/master)
[![coverage report](https://gitlab.com/agravgaard/cbctrecon/badges/master/coverage.svg)](https://gitlab.com/agravgaard/cbctrecon/commits/master)

**Proton dose calculation on scatter-corrected CBCT image: Feasibility study for adaptive proton therapy**
http://dx.doi.org/10.1118/1.4923179

**Supported OS/compilers**
 - Windows, 8.1 & 10
   - MSVC 2019
   - icc 2017 & 2019
 - Linux, Arch/Manjaro and Ubuntu
   - GCC 7, 8, 9
   - Clang 6, 9
   - icc 2019
 - MacOS, High Sierra
   - xcode

Some compilers and OSs are not tested as often.
I expect the code to compile on most systems with a C++14 compatible compiler.

In order to compile the software, you must have installed the following prerequisites.

**Prerequisites**
 - Git
 - Cmake >3.11
 - Qt 5.X
   - On Windows Set `CMAKE_PREFIX_PATH="[Qt install dir]/5.[X.Y]/msvc20[ZZ]_64/"`
 - OpenCL 
   - Use `-DRTK_USE_OPENCL=ON` to get OpenCL versions of the CUDA algorithms (prefer CUDA on nVidia platforms)

**Optional, but recommended, dependencies**
 - FFTW - Add directory to PATH
 - CUDA >8.0
   - C++ compiler compatible with your version of CUDA. (See [this table](https://gist.github.com/ax3l/9489132))

**Below is downloaded and compiled automatically:**
 - [DCMTK - Latest](https://github.com/DCMTK/DCMTK)
 - [ITK - Latest](https://github.com/InsightSoftwareConsortium/ITK)
 - [RTK - Latest](https://github.com/SimonRit/RTK)
 - [Plastimatch - Latest](https://gitlab.com/plastimatch/plastimatch)

*Before installation consider:
goPMC binary files and a Visual Studio 2013 Debug Compiler is needed for goPMC support and the goPMC app must be compiled seperately.*

If you want to use CUDA without nVidia hardware, then use the PGI compiler (pgc, pgcc) with nvcc from the CUDA SDK and copy the cuda dll's from the pgi bin directory to your working directory.

**How to build**
I recommend that you look at the short bash script, [ci/ci.sh](https://gitlab.com/agravgaard/cbctrecon/blob/master/ci/ci.sh), for an example of how to configure and build the project.
