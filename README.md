# cbctrecon
CBCT Reconstruction toolkit for VARIAN
This project is a fork of the Cone Beam Reconstruction project be Yang-Kyun Park et al.

<b>Proton dose calculation on scatter-corrected CBCT image: Feasibility study for adaptive proton therapy </b>
http://dx.doi.org/10.1118/1.4923179

The aim of this fork is to make the recontruction toolkit working with Varian CBCTs, as the project mentioned above only have worked with Elekta CBCTs.

For the moment only Windows 64bit with a CUDA compatible GPU is supported.

In order to compile the software, you must have installed the following prerequisites* and compiled with Visual Studio 2013 when compilation is needed:

<ul>
  <li>Cmake </li>
  <li>(Perl, Tcl, OpenSSL) </li>
  <li>Qt 4.X </li>
  <li>gdcm </li>
  <li>FFTW </li>
  <li>CUDA (CUDA might have a working FFTW included) </li>
  <li>DCMTK (Turn DCMTK_OVERWRITE_WIN32_COMPILER_FLAGS off)</li>
  <li>VTK </li>
  <li>ITK </li>
  <li>RTK 1.0rc03 - modified**</li>
  <li>Plastimatch - modified***</li>
</ul>

<ul>
  <li>*   If version not explicitly stated, up-to-date version should work</li>
  <li>**  Additional and modified files will be added to the project folder..</li>
  <li>*** Additional and modified files will be added to the project folder..</li>
</ul>

andreasg@phys.au.dk
