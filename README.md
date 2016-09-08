# cbctrecon
CBCT Reconstruction toolkit for VARIAN
This project is a fork of the Cone Beam Reconstruction project be Yang-Kyun Park et al.

<b>Proton dose calculation on scatter-corrected CBCT image: Feasibility study for adaptive proton therapy </b>
http://dx.doi.org/10.1118/1.4923179

The aim of this fork is to make the recontruction toolkit working with Varian CBCTs, as the project mentioned above only have worked with Elekta CBCTs.
HNC compatibility added via Geoff Hugo's RTK fork.

For the moment Windows 64bit is supported.

In order to compile the software, you must have installed the following prerequisites* and compiled with Visual Studio 2013 when compilation is needed:

<ul>
  <li>Cmake </li>
  <li>Qt 5.X [Updated after recent commit]</li>
  <li>gdcm (optional)</li>
  <li>FFTW (optional but should be preferred in ITK. While cuFFT in RTK)</li>
  <li>CUDA [Only the cudart library is mandatory]</li>
  <li>DCMTK (Turn DCMTK_OVERWRITE_WIN32_COMPILER_FLAGS off)</li>
  <li>VTK - Latest</li>
  <li>ITK - Latest</li>
  <li>RTK -  https://github.com/agravgaard/RTK</li>
  <li>Plastimatch -  https://gitlab.com/agravgaard/plastimatch</li>
  <li>Bzlib2 - https://github.com/philr/bzip2-windows/releases </li>
</ul>

Avoid using HDF5 when compiling dependencies or deal with the compile linking problems consequences yourself.
As a rule of thumb: if RTK and Plastimatch compiles this will too.

andreasg@phys.au.dk
