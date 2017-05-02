## Superbuild for RTK
set(proj RTK)

# Set dependency list
set(${proj}_DEPENDENCIES ITKv4)

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  message(FATAL_ERROR "Enabling ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj} is not supported !")
endif()

# Sanity checks
if(DEFINED RTK_DIR AND NOT EXISTS ${RTK_DIR})
  unset(RTK_DIR CACHE)
  find_package(RTK REQUIRED)
endif()

if(NOT DEFINED RTK_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  set(EXTERNAL_PROJECT_OPTIONAL_ARGS)
  if(CbctRecon_BUILD_DICOM_SUPPORT)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      -DRTK_USE_SYSTEM_DCMTK:BOOL=${RTK_USE_SYSTEM_DCMTK}
      -DDCMTK_DIR:PATH=${DCMTK_DIR}
      )
  endif()

  if(NOT DEFINED git_protocol)
    set(git_protocol "git")
  endif()


  IF (USE_CUDA)
    SET(RTK_USE_CUDA ON)
  ELSE (USE_CUDA)
    SET(RTK_USE_CUDA OFF)
  ENDIF (USE_CUDA)

  IF (USE_OPENCL)
    SET(RTK_USE_OPENCL ON)
    SET(RTK_REPO "agravgaard")
  ELSE (USE_OPENCL)
    SET(RTK_USE_OPENCL OFF)
    SET(RTK_REPO "SimonRit")
  ENDIF (USE_OPENCL)

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY
    "${git_protocol}://github.com/${RTK_REPO}/RTK.git"
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
      -DRTK_INSTALL_BIN_DIR:PATH=${CbctRecon_INSTALL_BIN_DIR}
      -DRTK_INSTALL_LIB_DIR:PATH=${CbctRecon_INSTALL_LIB_DIR}
      -DRTK_USE_GIT_PROTOCOL:BOOL=${CbctRecon_USE_GIT_PROTOCOL}
      -DRTK_USE_SYSTEM_ITK:BOOL=${RTK_USE_SYSTEM_ITK}
      -DITK_DIR:PATH=${ITK_DIR}
      -DCUDA_HAVE_GPU:BOOL=${RTK_USE_CUDA}
      -DRTK_USE_CUDA:BOOL=${RTK_USE_CUDA}
      -DOPENCL_HAVE_GPU:BOOL=${RTK_USE_OPENCL} # OpenCL flags is ignored in SimonRit version.
      -DRTK_USE_OPENCL:BOOL=${RTK_USE_OPENCL}
      -DOPENCL_INCLUDE_DIRS:PATH=${OpenCL_INCLUDE_DIR}
      -DOPENCL_LIBRARIES:FILEPATH=${OpenCL_LIBRARY}
      -DUSE_BZIP2:BOOL=${USE_BZIP2}
      -DGIT_EXECUTABLE:FILEPATH=${GIT_EXECUTABLE}
      ${EXTERNAL_PROJECT_OPTIONAL_ARGS}
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )

  ExternalProject_GenerateProjectDescription_Step(${proj})

  set(RTK_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  #-----------------------------------------------------------------------------
  # Launcher setting specific to build tree

  set(_lib_subdir lib)
  if(WIN32)
    set(_lib_subdir bin)
  endif()

  # library paths
  set(${proj}_LIBRARY_PATHS_LAUNCHER_BUILD ${RTK_DIR}/${_lib_subdir}/<CMAKE_CFG_INTDIR>)
  mark_as_superbuild(
    VARS ${proj}_LIBRARY_PATHS_LAUNCHER_BUILD
    LABELS "LIBRARY_PATHS_LAUNCHER_BUILD"
    )

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

mark_as_superbuild(
  VARS RTK_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
