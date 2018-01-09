# cbctrecon
CBCT Reconstruction toolkit for Elekta and Varian type cone beam projections

**Proton dose calculation on scatter-corrected CBCT image: Feasibility study for adaptive proton therapy**
http://dx.doi.org/10.1118/1.4923179

For the moment Windows 64bit is supported. *Linux (Arch/Manjaro) support is under development*
MacOS build support has just been added (CUDA not testet)

In order to compile the software, you must have installed the following prerequisites and a C++ compiler compatible with your version of CUDA.

**Prerequisites**
 - Git</li>
 - Cmake - Latest
 - Qt 5.X - On Windows Set `Qt5_DIR="[Qt install dir]/5.[X.Y]/msvc20[ZZ]\_64/lib/cmake/Qt5"`
 - FFTW (optional but recommended) - Add directory to PATH and copy dll's after compilation
 - CUDA and/or an OpenCL SDK

**Below is downloaded and compiled automatically:**
 - DCMTK - Latest: git://dcmtk.org/dcmtk
 - ITK - Latest: https://github.com/InsightSoftwareConsortium/ITK
 - RTK - Latest: https://github.com/SimonRit/RTK
 - Plastimatch - https://gitlab.com/plastimatch/plastimatch


*Before installation consider:
goPMC binary files and a Visual Studio 2013 Debug Compiler is needed for goPMC support and the goPMC app must be compiled seperately.*

If you want to use Visual Studio 2015, you must use CUDA 8.0 or above due to compiler incompatibilities.

If you want to use Visual Studio 2017, you must use CUDA 9.0 or above AND set `CUDA_HOST_COMPILER="C:/Program Files (x86)/Microsoft Visual Studio/[VS edition]/VC/Tools/MSVC/[cl version]/bin/Hostx64/x64/cl.exe"` due to compiler incompatibilities

If you want to use CUDA without nVidia hardware, then use the PGI compiler (pgc, pgcc) with nvcc from the CUDA SDK and copy the cuda dll's from the pgi bin directory to your working directory.

If you want to use the PGI compiler, remember to match the cuda versions i.e. update to version 16 to use CUDA 8.0

Avoid using HDF5 when compiling dependencies or deal with the compile linking problems yourself.

For older versions of Qt, ITK, CUDA check the RTK and Plastimatch compilation instructions (Good luck).
