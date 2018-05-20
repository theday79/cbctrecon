SET(DEPENDENCIES DCMTK ITKv4 Plastimatch)
IF (NOT ${RTK_AS_ITK_EXTERNAL})
list(APPEND DEPENDENCIES RTK)
ENDIF()

mark_as_superbuild(DEPENDENCIES:STRING)
include(ExternalProjectAddSource)

ExternalProject_Include_Dependencies(CbctRecon DEPENDS_VAR DEPENDENCIES)


#------------------------------------------------------------------------------
# Configure and build CbctRecon
#------------------------------------------------------------------------------
set(proj CbctRecon)

ExternalProject_Add(${proj}
  ${${proj}_EP_ARGS}
  DEPENDS ${DEPENDENCIES}
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR ${CbctRecon_BINARY_INNER_SUBDIR}
  DOWNLOAD_COMMAND ""
  UPDATE_COMMAND ""
  CMAKE_CACHE_ARGS
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
	  -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
    -DADDITIONAL_C_FLAGS:STRING=${ADDITIONAL_C_FLAGS}
    -DADDITIONAL_CXX_FLAGS:STRING=${ADDITIONAL_CXX_FLAGS}
    -DCbctRecon_REQUIRED_C_FLAGS:STRING=${CbctRecon_REQUIRED_C_FLAGS}
    -DCbctRecon_REQUIRED_CXX_FLAGS:STRING=${CbctRecon_REQUIRED_CXX_FLAGS}
    -DCbctRecon_SUPERBUILD:BOOL=OFF
    -DCbctRecon_SUPERBUILD_DIR:PATH=${CbctRecon_BINARY_DIR}
    -D${CbctRecon_MAIN_PROJECT}_APPLICATION_NAME:STRING=${${CbctRecon_MAIN_PROJECT}_APPLICATION_NAME}
    -D${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_MAJOR:STRING=${${CbctReconAppName}_VERSION_MAJOR}
    -D${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_MINOR:STRING=${${CbctReconAppName}_VERSION_MINOR}
    -D${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_PATCH:STRING=${${CbctReconAppName}_VERSION_PATCH}
    -D${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_TWEAK:STRING=${${CbctReconAppName}_VERSION_TWEAK}
    -D${CbctRecon_MAIN_PROJECT_APPLICATION_NAME}_VERSION_RC:STRING=${${CbctReconAppName}_VERSION_RC}
    -DUSE_CUDA:BOOL=${USE_CUDA}
    -DCUDA_SEPARABLE_COMPILATION:BOOL=${CUDA_SEPARABLE_COMPILATION}
    -DUSE_OPENCL:BOOL=${USE_OPENCL}
    -DUSE_OPENMP:BOOL=${USE_OPENMP}
    -DUSE_CLFFT:BOOL=${USE_CLFFT}
    -DUSE_GPMC:BOOL=${USE_GPMC}
	  -DUSE_LOWPASS_FFT:BOOL=${USE_LOWPASS_FFT}
    -DFFTW_FOUND:BOOL=${FFTW_FOUND}
    -DFFTW_DIR:PATH=${FFTW_DIR}
    -DOpenCL_INCLUDE_DIR:PATH=${OpenCL_INCLUDE_DIR}
    -DOpenCL_LIBRARY:FILEPATH=${OpenCL_LIBRARY}
    -DCLFFT_FOUND:BOOL=${CLFFT_FOUND}
    -D_CLFFT_LIBRARY:PATH=${_CLFFT_LIBRARY}
    -D_CLFFT_INCLUDE_DIRS:PATH=${_CLFFT_INCLUDE_DIRS}
    -DCLFFT_ROOT_DIR:PATH=${CLFFT_ROOT_DIR}
    -DQt5_DIR:PATH=${Qt5_DIR}
    -DITK_DIR:PATH=${ITK_DIR}
    -DRTK_AS_ITK_EXTERNAL:BOOL=${RTK_AS_ITK_EXTERNAL}
    -DRTK_DIR:PATH=${RTK_DIR}
    -DDCMTK_DIR:PATH=${DCMTK_DIR}
    -DPlastimatch_DIR:PATH=${Plastimatch_DIR}
    -DPlastimatch_BUILD_DIR:PATH=${Plastimatch_BUILD_DIR}
    -DPlastimatch_SRC:PATH=${Plastimatch_SRC}
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
    -DCbctRecon_APPLICATIONS_DIR:PATH=${CbctRecon_APPLICATIONS_DIR}
    -DCbctRecon_EXTENSION_SOURCE_DIRS:STRING=${CbctRecon_EXTENSION_SOURCE_DIRS}
    -DExternalData_OBJECT_STORES:PATH=${ExternalData_OBJECT_STORES}
    ${EXTERNAL_PROJECT_OPTIONAL_ARGS}
  )

# This custom external project step forces the build and later
# steps to run whenever a top level build is done...
#
# BUILD_ALWAYS flag is available in CMake 3.1 that allows force build
# of external projects without this workaround. Remove this workaround
# and use the CMake flag instead, when CbctRecon's required minimum CMake
# version will be at least 3.1.
#
if(CMAKE_CONFIGURATION_TYPES)
  set(BUILD_STAMP_FILE "${CMAKE_CURRENT_BINARY_DIR}/${proj}-prefix/src/${proj}-stamp/${CMAKE_CFG_INTDIR}/${proj}-build")
else()
  set(BUILD_STAMP_FILE "${CMAKE_CURRENT_BINARY_DIR}/${proj}-prefix/src/${proj}-stamp/${proj}-build")
endif()
ExternalProject_Add_Step(${proj} forcebuild
  COMMAND ${CMAKE_COMMAND} -E remove ${BUILD_STAMP_FILE}
  COMMENT "Forcing build step for '${proj}'"
  DEPENDEES build
  ALWAYS 1
  )
