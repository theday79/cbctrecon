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

  IF (USE_OPENCL)
    SET(Plastimatch_DISABLE_OPENCL OFF)
  ELSE (USE_OPENCL)
    SET(Plastimatch_DISABLE_OPENCL ON)
  ENDIF (USE_OPENCL)

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY
    "${git_protocol}://github.com/agravgaard/plastimatch.git"
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
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DADDITIONAL_C_FLAGS:STRING=${ADDITIONAL_C_FLAGS}
      -DADDITIONAL_CXX_FLAGS:STRING=${ADDITIONAL_CXX_FLAGS}
      -DBUILD_TESTING:BOOL=OFF
      -DPLM_INSTALL_BIN_DIR:STRING=${CbctRecon_INSTALL_BIN_DIR}
      -DPLM_INSTALL_LIB_DIR:STRING=${CbctRecon_INSTALL_LIB_DIR}
      -DPLM_USE_GIT_PROTOCOL:BOOL=${CbctRecon_USE_GIT_PROTOCOL}
      -DDCMTK_DIR:STRING=${DCMTK_DIR}
      -DITK_DIR:STRING=${ITK_DIR}
      -DZLIB_DIR:STRING=${ZLIB_DIR}
      -DPLM_CONFIG_DISABLE_CUDA:BOOL=${Plastimatch_DISABLE_CUDA}
      -DPLM_CONFIG_DISABLE_OPENCL:BOOL=${Plastimatch_DISABLE_OPENCL}
      -DPLM_CONFIG_DISABLE_REG23:BOOL=ON
      -DOPENCL_INCLUDE_DIRS:STRING=${OpenCL_INCLUDE_DIR}
      -DOPENCL_LIBRARIES:FILEPATH=${OpenCL_LIBRARY}
      -DGIT_EXECUTABLE:FILEPATH=${GIT_EXECUTABLE}
      ${EXTERNAL_PROJECT_OPTIONAL_ARGS}
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )

  ExternalProject_GenerateProjectDescription_Step(${proj})

  set(Plastimatch_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  #-----------------------------------------------------------------------------
  # Launcher setting specific to build tree

  # library paths
  set(${proj}_LIBRARY_PATHS_LAUNCHER_BUILD ${Plastimatch_DIR}/Plastimatch-build/bin/<CMAKE_CFG_INTDIR>)
  if(CbctRecon_USE_QtTesting)
    list(APPEND ${proj}_LIBRARY_PATHS_LAUNCHER_BUILD
      ${Plastimatch_DIR}/QtTesting-build/<CMAKE_CFG_INTDIR>
      )
  endif()
  mark_as_superbuild(
    VARS ${proj}_LIBRARY_PATHS_LAUNCHER_BUILD
    LABELS "LIBRARY_PATHS_LAUNCHER_BUILD"
    )

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

set(Plastimatch_BUILD_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
set(Plastimatch_SRC ${CMAKE_BINARY_DIR}/${proj})

#mark_as_superbuild(
#  VARS Plastimatch_DIR:PATH
#  LABELS "FIND_PACKAGE"
#  )
