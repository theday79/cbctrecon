OPTION(USE_GPMC "Enable gpu Proton Monte Carlo (gPMC) dose re-calculation" ON)
IF(USE_GPMC)
  SET(MEMORY_SIZE 32 CACHE STRING "Computer memory in GB (for ITK)")
  include(ExternalProject)
  ExternalProject_Add(gPMC ## Yes, external is needed because otherwise dlls will interfere with local DCMTK
	DOWNLOAD_COMMAND ""    ## Also, DCMTK must be installed in a path not called %PROGRAMFILES%/DCMTK !!!
	SOURCE_DIR "${PROJECT_SOURCE_DIR}/gPMC"
	# CONFIGURE_COMMAND ""
	CMAKE_GENERATOR "Visual Studio 12 2013 Win64"
	CMAKE_GENERATOR_TOOLSET "v120" ## Visual Studio 2013
	CONFIGURATIONS "Debug" ## Release has not worked yet, probably some "if debug" directives in goPMC or dependencies
	INSTALL_DIR "${CMAKE_BINARY_DIR}/install"
	CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install
      -DMEMORY_SIZE:STRING=${MEMORY_SIZE}
  )
  ADD_DEFINITIONS(-DUSE_GPMC=TRUE)
ENDIF(USE_GPMC)
