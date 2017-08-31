goPMC

This package consists of interfaces and dynamic link libraries (dll) of a Monte Carlo 
simulation engine for proton therapy, goPMC. The binaries were built with Visual Studio 
2013 in debug x64 configuration. Undefined behaviors may occur when different compilers 
or configurations are used.

System requirements:
OpenCL 1.1 or higher.
CPU/GPU supporting OpenCL.

Subdirectories:
include: 
	goPMC.h: goPMC interface.
	cl.hpp: c++ wrapper for OpenCL 1.1.
bin:
	dcmtk.dll(.lib) and libDicomRT.dll(.lib): Libraries for reading and processing Dicom CT data.
	goPMC.dll(.lib): goPMC library.
input: 
	Physics input data.




Related publications:

[1] GPU-based fast Monte Carlo dose calculation for proton therapy
Xun Jia, Jan Schuemann, Harald Paganetti and Steve B. Jiang
Physics in Medicine and Biology, Volume 57, Number 23

[2] Validation of a GPU-based Monte Carlo code (gPMC) for proton radiation therapy: clinical cases study
Drosoula Giantsoudi, Jan Schuemann, Xun Jia, Stephen Dowdell, Steve B. Jiang and Harald Paganetti
Physics in Medicine and Biology, Volume 60, Number 6

[3] Recent developments and comprehensive evaluations of a GPU-based Monte Carlo package for proton therapy
Nan Qin, Pablo Botas, Drosoula Giantsoudi, Jan Schuemann, Zhen Tian, Steve B. Jiang, Harald Paganetti and Xun Jia
Accepted by Physics in Medicine and Biology


For more questions please email to nan.qin@utsouthwestern.edu.