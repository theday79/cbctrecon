OPTION(USE_GPMC "Enable gpu Proton Monte Carlo (gPMC) dose re-calculation" ON)
IF(USE_GPMC)
  SET(proj gPMC)

  # Set dependency list
  set(${proj}_DEPENDENCIES ITKv4)

  # Include dependent projects if any
  ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

  SET(MEMORY_SIZE 32 CACHE STRING "Computer memory in GB (for ITK)")
  ExternalProject_Add(${proj} ## Yes, external is needed because otherwise dlls will interfere with local DCMTK
    DEPENDS ${${proj}_DEPENDENCIES}            ## Because download should only happen once.
	  DOWNLOAD_COMMAND ""      ## Also, DCMTK must be installed in a path not called %PROGRAMFILES%/DCMTK !!!
	  SOURCE_DIR "${PROJECT_SOURCE_DIR}/${proj}"
	  # CONFIGURE_COMMAND ""
	  CMAKE_GENERATOR "Visual Studio 12 2013 Win64"
	  CMAKE_GENERATOR_TOOLSET "v120" ## Visual Studio 2013
	  CONFIGURATIONS "Debug" ## Release has not worked yet, probably some "if debug" directives in goPMC or dependencies
	  INSTALL_DIR "${CMAKE_BINARY_DIR}/install"
	  CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install
      -DMEMORY_SIZE:STRING=${MEMORY_SIZE}
      -DITK_DIR=${ITK_DIR}
  )
  ADD_DEFINITIONS(-DUSE_GPMC=TRUE)

  ExternalProject_GenerateProjectDescription_Step(${proj})

  set(${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  #-----------------------------------------------------------------------------
  # Launcher setting specific to build tree

  mark_as_superbuild(
    VARS gPMC_DIR:PATH
    LABELS "FIND_PACKAGE"
    )

ENDIF(USE_GPMC)
