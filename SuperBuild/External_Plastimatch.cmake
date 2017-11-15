## Superbuild for Plastimatch
set(proj Plastimatch)

# Set dependency list
set(${proj}_DEPENDENCIES ITKv4 DCMTK)

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  message(FATAL_ERROR "Enabling ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj} is not supported !")
endif()

# Sanity checks
if(DEFINED Plastimatch_DIR AND NOT EXISTS ${Plastimatch_DIR})
  unset(Plastimatch_DIR CACHE)
  find_package(Plastimatch REQUIRED)
endif()

if(NOT DEFINED Plastimatch_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  set(EXTERNAL_PROJECT_OPTIONAL_ARGS)
  if(CbctRecon_BUILD_DICOM_SUPPORT)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      -DPlastimatch_USE_SYSTEM_DCMTK:BOOL=${Plastimatch_USE_SYSTEM_DCMTK}
      -DDCMTK_DIR:PATH=${DCMTK_DIR}
      )
  endif()

  if(NOT DEFINED git_protocol)
    set(git_protocol "git")
  endif()


  IF (USE_CUDA)
    SET(Plastimatch_DISABLE_CUDA OFF)
  ELSE (USE_CUDA)
    SET(Plastimatch_DISABLE_CUDA ON)
  ENDIF (USE_CUDA)

  IF (OpenCL_FOUND)
    SET(Plastimatch_DISABLE_OPENCL OFF)
  ELSE (OpenCL_FOUND)
    SET(Plastimatch_DISABLE_OPENCL ON)
  ENDIF (OpenCL_FOUND)

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY
    #"${git_protocol}://github.com/agravgaard/plastimatch.git"
    "https://gitlab.com/plastimatch/plastimatch.git"
    QUIET
    )

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG
    "b8587e55ee7d6843a6ab2688432ade35fc7581b6"
    QUIET
    )

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY "${${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY}"
    # GIT_TAG "${${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG}"
    SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj}
    BINARY_DIR ${proj}-build
    CMAKE_CACHE_ARGS
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      -DCMAKE_ECLIPSE_VERSION:STRING=${CMAKE_ECLIPSE_VERSION}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS:STRING=${PLAST_CXX_FLAGS}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DCMAKE_LINKER:FILEPATH=${CMAKE_LINKER}
      -DADDITIONAL_C_FLAGS:STRING=${ADDITIONAL_C_FLAGS}
      -DADDITIONAL_CXX_FLAGS:STRING=${ADDITIONAL_CXX_FLAGS}
      -DBUILD_TESTING:BOOL=${BUILD_TESTING}
      -DPLM_INSTALL_BIN_DIR:PATH=${CbctRecon_INSTALL_BIN_DIR}
      -DPLM_INSTALL_LIB_DIR:PATH=${CbctRecon_INSTALL_LIB_DIR}
      -DPLM_USE_GIT_PROTOCOL:BOOL=${CbctRecon_USE_GIT_PROTOCOL}
      -DDCMTK_DIR:PATH=${DCMTK_DIR}
      -DITK_DIR:PATH=${ITK_DIR}
      -DZLIB_DIR:PATH=${ZLIB_DIR}
      -DPLM_CONFIG_DISABLE_CUDA:BOOL=${Plastimatch_DISABLE_CUDA}
      -DPLM_CONFIG_DISABLE_OPENCL:BOOL=${Plastimatch_DISABLE_OPENCL}
      -DPLM_CONFIG_DISABLE_REG23:BOOL=ON
      -DOPENCL_INCLUDE_DIRS:PATH=${OpenCL_INCLUDE_DIR}
      -DOPENCL_LIBRARIES:FILEPATH=${OpenCL_LIBRARY}
	  -DCUDA_TOOLKIT_ROOT_DIR:PATH=${CUDA_TOOLKIT_ROOT_DIR}
      -DFFTW_FOUND:BOOL=${FFTW_FOUND}
      -DFFTW_LIBRARIES:PATH=${FFTW_LIBRARIES}
      -DFFTW_INCLUDE_DIR:PATH=${FFTW_INCLUDE_DIR}
      -DFFTW_DIR:PATH=${FFTW_INCLUDE_DIR}
      -DGIT_EXECUTABLE:FILEPATH=${GIT_EXECUTABLE}
      # TBB
      -DTBB_INCLUDE_DIR:PATH=${TBB_INCLUDE_DIR}
      -DTBB_LIBRARY_DEBUG:PATH=${TBB_LIBRARY_DEBUG}
      -DTBB_LIBRARY_RELEASE:PATH=${TBB_LIBRARY_RELEASE}
      -DTBB_MALLOC_INCLUDE_DIR:PATH=${TBB_MALLOC_INCLUDE_DIR}
      -DTBB_MALLOC_LIBRARY_DEBUG:PATH=${TBB_MALLOC_LIBRARY_DEBUG}
      -DTBB_MALLOC_LIBRARY_RELEASE:PATH=${TBB_MALLOC_LIBRARY_RELEASE}
      -DTBB_MALLOC_PROXY_INCLUDE_DIR:PATH=${TBB_MALLOC_PROXY_INCLUDE_DIR}
      -DTBB_MALLOC_PROXY_LIBRARY_DEBUG:PATH=${TBB_MALLOC_PROXY_LIBRARY_DEBUG}
      -DTBB_MALLOC_PROXY_LIBRARY_RELEASE:PATH=${TBB_MALLOC_PROXY_LIBRARY_RELEASE}
      ${EXTERNAL_PROJECT_OPTIONAL_ARGS}
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/${proj}-install
    # INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )

  ExternalProject_GenerateProjectDescription_Step(${proj})

  set(Plastimatch_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
  set(Plastimatch_BUILD_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
  set(Plastimatch_SRC ${CMAKE_BINARY_DIR}/${proj})

  #-----------------------------------------------------------------------------
  # Launcher setting specific to build tree

  set(_lib_subdir lib)
  if(WIN32)
    set(_lib_subdir bin)
  endif()

  # library paths
  set(${proj}_LIBRARY_PATHS_LAUNCHER_BUILD ${Plastimatch_DIR}/${_lib_subdir}/<CMAKE_CFG_INTDIR>)
  mark_as_superbuild(
    VARS ${proj}_LIBRARY_PATHS_LAUNCHER_BUILD
    LABELS "LIBRARY_PATHS_LAUNCHER_BUILD"
    )

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

if(${TBB_FOUND})
  add_custom_command(
    TARGET ${proj} PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${TBB_LIBRARY_RELEASE}
    ${plastimatch_DIR}/${CMAKE_CFG_INTDIR}
  )
endif()

mark_as_superbuild(
  VARS Plastimatch_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
