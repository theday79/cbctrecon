# cbctrecon
CBCT Reconstruction toolkit for Elekta and Varian type cone beam projections

<b>Proton dose calculation on scatter-corrected CBCT image: Feasibility study for adaptive proton therapy </b>
http://dx.doi.org/10.1118/1.4923179

For the moment Windows 64bit is supported
In order to compile the software, you must have installed the following prerequisites and compiled with Visual Studio 2013 (or 2015 with CUDA 8.0 or above) when compilation is needed:

<ul>
  <li>Git</li>
  <li>Cmake - Latest</li>
  <li>Qt 5.X - On Windows Set Qt5_DIR="[Qt install dir]/5.[X.Y]/msvc20[ZZ]_64/lib/cmake/Qt5"</li>
  <li>FFTW (optional but recommended) - Add directory to PATH and copy dll's after compilation</li>
  <li>CUDA and/or an OpenCL SDK</li>
  <li><b>Below is downloaded and compiled automatically:</b></li>
  <li>DCMTK - Latest: git://dcmtk.org/dcmtk </li>
  <li>ITK - Latest: https://github.com/InsightSoftwareConsortium/ITK </li>
  <li>RTK - Latest: https://github.com/SimonRit/RTK </li>
  <li>Plastimatch - https://gitlab.com/plastimatch/plastimatch </li>
</ul>

*Before installation consider:
goPMC binary files and a Visual Studio 2013 Debug Compiler is needed for goPMC support and the goPMC app must be compiled seperately.

If you want to use Visual Studio 2015, you must use CUDA 8.0 or above due to compiler incompatibilities
If you want to use Visual Studio 2017, you must use CUDA 9.0 or above AND set CUDA_HOST_COMPILER="C:/Program Files (x86)/Microsoft Visual Studio/[VS edition]/VC/Tools/MSVC/[cl version]/bin/Hostx64/x64/cl.exe" due to compiler incompatibilities

If you want to use CUDA without nVidia hardware, then use the PGI compiler (pgc, pgcc) with nvcc from the CUDA SDK and copy the cuda dll's from the pgi bin directory to your working directory.
If you want to use the PGI compiler, remember to match the cuda versions i.e. update to version 16 to use CUDA 8.0
Avoid using HDF5 when compiling dependencies or deal with the compile linking problems yourself.
For older versions of Qt, ITK, CUDA check the RTK and Plastimatch compilation instructions (Good luck).
