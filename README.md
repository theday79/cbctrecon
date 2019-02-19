# CbctRecon
CBCT Reconstruction toolkit for Elekta and Varian type cone beam projections

[![pipeline status](https://gitlab.com/agravgaard/cbctrecon/badges/master/pipeline.svg)](https://gitlab.com/agravgaard/cbctrecon/commits/WEPLvisualization)
[![Build status](https://ci.appveyor.com/api/projects/status/ek8g59aguufokj3q?svg=true)](https://ci.appveyor.com/project/agravgaard/cbctrecon)
[![Coverage Status](https://coveralls.io/repos/gitlab/agravgaard/cbctrecon/badge.svg?branch=master)](https://coveralls.io/gitlab/agravgaard/cbctrecon?branch=master)

**Proton dose calculation on scatter-corrected CBCT image: Feasibility study for adaptive proton therapy**
http://dx.doi.org/10.1118/1.4923179

Software will be tested on Linux (Arch, Manjaro and Ubuntu), MacOS (High Sierra), Windows (8.1 , 10 and Server 2012).
Using GCC (7, 8), Clang (6, 7 and latest xcode), ICC (17, 18), MSVC (2015 and 2017).
((Currently Clang is unable to compile anything requiring pthread.h or nanosleep on my system and will not be tested until it is fixed))

Some compilers and OSs are not tested as often.

In order to compile the software, you must have installed the following prerequisites and a C++ compiler compatible with your version of CUDA. (See [this table](https://gist.github.com/ax3l/9489132))

**Prerequisites**
 - Git
 - Cmake - Latest
 - Qt 5.X - On Windows Set `Qt5_DIR="[Qt install dir]/5.[X.Y]/msvc20[ZZ]_64/lib/cmake/Qt5"`
 - FFTW (optional but recommended) - Add directory to PATH
 - CUDA and/or an OpenCL SDK

**Below is downloaded and compiled automatically:**
 - DCMTK - Latest: https://github.com/commontk/DCMTK
 - ITK - 4.13.1: https://github.com/InsightSoftwareConsortium/ITK
 - RTK - Latest: https://github.com/SimonRit/RTK
 - Plastimatch - latest - https://gitlab.com/plastimatch/plastimatch

**Optional dependencies (downloaded and compiled automatically):**
 - tinyrefl - Latest: https://gitlab.com/Manu343726/tinyrefl

*Before installation consider:
goPMC binary files and a Visual Studio 2013 Debug Compiler is needed for goPMC support and the goPMC app must be compiled seperately.*

Plastimatch does not yet support ITK 5, but we apply a patch as a workaround until official support.

If you want to use Visual Studio 2015, you must use CUDA 8.0 or above due to compiler incompatibilities.

If you want to use Visual Studio 2017 (v14.11 or **earlier**), you must use CUDA 9.0 or above AND set `CUDA_HOST_COMPILER="C:/Program Files (x86)/Microsoft Visual Studio/[VS edition]/VC/Tools/MSVC/[cl version]/bin/Hostx64/x64/cl.exe"` due to compiler incompatibilities

If you want to use CUDA without nVidia hardware, then use the PGI compiler (pgc, pgcc) with nvcc from the CUDA SDK and copy the cuda dll's from the pgi bin directory to your working directory.

Avoid using Intel TBB when compiling dependencies or deal with the compile linking problems yourself.

For older versions of Qt, ITK, CUDA check the RTK and Plastimatch compilation instructions (Good luck).
