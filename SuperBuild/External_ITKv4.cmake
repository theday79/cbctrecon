
set(proj ITKv4)

# Set dependency list
set(${proj}_DEPENDENCIES "DCMTK") # "zlib")

if(CbctRecon_BUILD_ITKPython)
  list(APPEND ${proj}_DEPENDENCIES Swig python)
endif()

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  unset(ITK_DIR CACHE)
  find_package(ITK REQUIRED)
endif()

# Sanity checks
if(DEFINED ITK_DIR AND NOT EXISTS ${ITK_DIR})
  message(FATAL_ERROR "ITK_DIR variable is defined but corresponds to nonexistent directory")
endif()

if(NOT DEFINED ITK_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  if(NOT DEFINED git_protocol)
      set(git_protocol "git")
  endif()

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY
    "${git_protocol}://github.com/InsightSoftwareConsortium/ITK.git"
    QUIET
    )

  # ITK release v4.11.0 from 2017.01.22 with
  # * CbctRecon patches for CMP0042
  # * Revert "BUG: Binary include directory added only if enabled shared" (ITK-3529)
  # * Revert "BUG: Include installed modules in ITK_MODULES_ENABLED" (ITK-3529)
  # * COMP: Add missing headers in OpenCVImageBridge
  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG
    "d92873e33e8a54e933e445b92151191f02feab42" #"v4.13.0"
    # "1619816e48b327c4c486b76aef0c229af65b92b2" # CbctRecon-v4.11.0-2017-01-22
    QUIET
    )

  set(EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS)

  if(CbctRecon_USE_PYTHONQT OR CbctRecon_BUILD_ITKPython)
    # XXX Ensure python executable used for ITKModuleHeaderTest
    #     is the same as CbctRecon.
    #     This will keep the sanity check implemented in CbctReconConfig.cmake
    #     quiet.
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
      -DPYTHON_EXECUTABLE:PATH=${PYTHON_EXECUTABLE}
      )
  endif()

  if(CbctRecon_BUILD_ITKPython)

    # Custom name for the components associated with ITK
    # wrapping install rules enabling CbctRecon to optionally
    # package ITK Wrapping in CbctRecon installer by simply
    # toggling the CbctRecon_INSTALL_ITKPython option.
    set(CbctRecon_WRAP_ITK_INSTALL_COMPONENT_IDENTIFIER "Wrapping")
    mark_as_superbuild(CbctRecon_WRAP_ITK_INSTALL_COMPONENT_IDENTIFIER:STRING)

    set(PY_SITE_PACKAGES_PATH lib/Python/${pythonpath_subdir}/site-packages)

    list(APPEND EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
      -DPYTHON_LIBRARY:FILEPATH=${PYTHON_LIBRARY}
      -DPYTHON_INCLUDE_DIR:PATH=${PYTHON_INCLUDE_DIR}
      -DSWIG_EXECUTABLE:PATH=${SWIG_EXECUTABLE}
      -DITK_USE_SYSTEM_SWIG:BOOL=ON
      -DITK_LEGACY_SILENT:BOOL=ON
      -DWRAP_ITK_INSTALL_COMPONENT_IDENTIFIER:STRING=${CbctRecon_WRAP_ITK_INSTALL_COMPONENT_IDENTIFIER}
      -DPY_SITE_PACKAGES_PATH:STRING=${PY_SITE_PACKAGES_PATH}
      )
  endif()


  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY "${${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY}"
    GIT_TAG "${${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG}"
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    CMAKE_CACHE_ARGS
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      -DCMAKE_ECLIPSE_VERSION:STRING=${CMAKE_ECLIPSE_VERSION}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
	  -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
      -DCMAKE_LINKER:FILEPATH=${CMAKE_LINKER}
      -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
      -DITK_INSTALL_ARCHIVE_DIR:PATH=${CbctRecon_INSTALL_LIB_DIR}
      -DITK_INSTALL_LIBRARY_DIR:PATH=${CbctRecon_INSTALL_LIB_DIR}
      -DBUILD_TESTING:BOOL=${BUILD_TESTING}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DITK_BUILD_DEFAULT_MODULES:BOOL=ON
      -DModule_ITKReview:BOOL=ON
      -DITK_WRAPPING:BOOL=OFF
      -DITK_WRAP_PYTHON:BOOL=${CbctRecon_BUILD_ITKPython}
      -DExternalData_OBJECT_STORES:PATH=${ExternalData_OBJECT_STORES}
      # VTK
      -DModule_ITKVtkGlue:BOOL=OFF
      # DCMTK
      -DITK_USE_SYSTEM_DCMTK:BOOL=ON
      -DDCMTK_DIR:PATH=${DCMTK_DIR}
      -DModule_ITKIODCMTK:BOOL=${CbctRecon_BUILD_DICOM_SUPPORT}
      # ZLIB
      -DITK_USE_SYSTEM_ZLIB:BOOL=OFF
      # -DZLIB_ROOT:PATH=${ZLIB_ROOT}
      # -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
      # -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
      # TBB
      -DModule_TBBImageToImageFilter:BOOL=OFF #${TBB_FOUND}
      #-DTBB_INCLUDE_DIR:PATH=${TBB_INCLUDE_DIR}
      #-DTBB_LIBRARY_DEBUG:PATH=${TBB_LIBRARY_DEBUG}
      #-DTBB_LIBRARY_RELEASE:PATH=${TBB_LIBRARY_RELEASE}
      #-DTBB_MALLOC_INCLUDE_DIR:PATH=${TBB_MALLOC_INCLUDE_DIR}
      #-DTBB_MALLOC_LIBRARY_DEBUG:PATH=${TBB_MALLOC_LIBRARY_DEBUG}
      #-DTBB_MALLOC_LIBRARY_RELEASE:PATH=${TBB_MALLOC_LIBRARY_RELEASE}
      #-DTBB_MALLOC_PROXY_INCLUDE_DIR:PATH=${TBB_MALLOC_PROXY_INCLUDE_DIR}
      #-DTBB_MALLOC_PROXY_LIBRARY_DEBUG:PATH=${TBB_MALLOC_PROXY_LIBRARY_DEBUG}
      #-DTBB_MALLOC_PROXY_LIBRARY_RELEASE:PATH=${TBB_MALLOC_PROXY_LIBRARY_RELEASE}
      # FFTW
      -DITK_USE_FFTWF:BOOL=${FFTW_FOUND}
      -DITK_USE_FFTWD:BOOL=${FFTW_FOUND}
      -DITK_USE_SYSTEM_FFTW:BOOL=${FFTW_FOUND}
      -DFFTW_FOUND:BOOL=${FFTW_FOUND}
      -DFFTW_LIBRARIES:PATH=${FFTW_LIBRARIES}
      -DFFTWD_LIB:PATH=${FFTW_LIBRARY_DOUBLE}
      -DFFTWF_LIB:PATH=${FFTW_LIBRARY_FLOAT}
      -DFFTWD_THREADS_LIB:PATH=${FFTW_THREADS_LIBRARY_DOUBLE}
      -DFFTWF_THREADS_LIB:PATH=${FFTW_THREADS_LIBRARY_FLOAT}
      -DFFTW_INCLUDE_DIR:PATH=${FFTW_INCLUDE_DIR}
      -DFFTW_INCLUDE_PATH:PATH=${FFTW_INCLUDE_DIR}
      # GPU
      -DITK_USE_GPU:BOOL=ON
      -DOPENCL_INCLUDE_DIRS:PATH=${OpenCL_INCLUDE_DIR}
      -DOPENCL_LIBRARIES:FILEPATH=${OpenCL_LIBRARY}
      -DITK_COMPUTER_MEMORY_SIZE:STRING=32
      ${EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS}
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )

  ExternalProject_GenerateProjectDescription_Step(${proj})

  set(ITK_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
  set(ITK_SRC_DIR ${CMAKE_BINARY_DIR}/${proj})
  #-----------------------------------------------------------------------------
  # Launcher setting specific to build tree

  set(_lib_subdir lib)
  if(WIN32)
    set(_lib_subdir bin)
  endif()

  # library paths
  set(${proj}_LIBRARY_PATHS_LAUNCHER_BUILD ${ITK_DIR}/${_lib_subdir}/<CMAKE_CFG_INTDIR>)
  mark_as_superbuild(
    VARS ${proj}_LIBRARY_PATHS_LAUNCHER_BUILD
    LABELS "LIBRARY_PATHS_LAUNCHER_BUILD"
    )

  if(CbctRecon_BUILD_ITKPython)
    # pythonpath
    set(${proj}_PYTHONPATH_LAUNCHER_BUILD
      ${ITK_DIR}/Wrapping/Generators/Python/<CMAKE_CFG_INTDIR>
      ${ITK_DIR}/lib/<CMAKE_CFG_INTDIR>
      ${ITK_DIR}/lib
      )
    mark_as_superbuild(
      VARS ${proj}_PYTHONPATH_LAUNCHER_BUILD
      LABELS "PYTHONPATH_LAUNCHER_BUILD"
      )
  endif()

  #-----------------------------------------------------------------------------
  # Launcher setting specific to install tree

  # Since ITK Wrapping is installed in the CbctRecon standard site-packages
  # location, there is no need to specify custom setting for the install
  # case.

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

if(${TBB_FOUND})
  add_custom_command(
    TARGET ${proj} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${TBB_LIBRARY_RELEASE}
    ${ITK_DIR}/lib/${CMAKE_CFG_INTDIR}
  )
endif()

mark_as_superbuild(
  VARS ITK_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
