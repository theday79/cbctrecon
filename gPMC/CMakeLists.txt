CMAKE_MINIMUM_REQUIRED (VERSION 3.8)

PROJECT(gPMC)

SET(GPMC_VERSION_MAJOR "2")
SET(GPMC_VERSION_MINOR "0")
SET(GPMC_VERSION_PATCH "0")

SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

FIND_PACKAGE(OpenCL REQUIRED)
FIND_PACKAGE(Gengetopt)

################################ ITK Superbuild ###########################
find_package(Git)

if(NOT GIT_FOUND)
  message(ERROR "Cannot find git. git is required for Superbuild")
endif()

option( USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)

set(git_protocol "git")
if(NOT USE_GIT_PROTOCOL)
  set(git_protocol "http")
endif()

include( ExternalProject )

if( NOT ITK_TAG )
  # ITK master 2017-03-06
  set( ITK_TAG "08dd93f8844169d92f5649b342ca5a167cfd6849" )
endif()
SET(MEMORY_SIZE 32 CACHE STRING "Cumputer memory in GB (for ITK)")
set(ep_common_args
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
  -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
)
ExternalProject_Add(ITK
  GIT_REPOSITORY "${git_protocol}://github.com/InsightSoftwareConsortium/ITK.git"
  GIT_TAG "${ITK_TAG}"
  SOURCE_DIR ITK
  BINARY_DIR ITK-build
  CMAKE_GENERATOR ${CMAKE_GENERATOR}
  CMAKE_GENERATOR_TOOLSET ${CMAKE_GENERATOR_TOOLSET}
  CMAKE_ARGS
    ${ep_common_args}
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_TESTING:BOOL=OFF
    -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF
    -DITKGroup_Core:BOOL=OFF
    -DITKGroup_IO:BOOL=OFF
    -DModule_ITKReview:BOOL=ON
    -DModule_ITKDCMTK:BOOL=OFF
    -DITK_LEGACY_SILENT:BOOL=ON
    -DITK_USE_KWSTYLE:BOOL=OFF
    -DITK_COMPUTER_MEMORY_SIZE:STRING=${MEMORY_SIZE}
    -DExternalData_OBJECT_STORES:STRING=${ExternalData_OBJECT_STORES}
  INSTALL_COMMAND ""
  ## LOG_BUILD 0
)

set(ITK_DIR ${CMAKE_BINARY_DIR}/ITK-build)
################################# ITK Superbuild ###########################
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

SET(GPMC_LIB_DIR "${PROJECT_SOURCE_DIR}/bin")
SET(GPMC_LIBRARIES
  "${GPMC_LIB_DIR}/goPMC.lib"
  "${GPMC_LIB_DIR}/dcmtk.lib"
  "${GPMC_LIB_DIR}/libDicomRT.lib"
)

WRAP_GGO(gPMC_GGO_C gPMC.ggo gPMC3Doutputdose_section.ggo gPMCVersion.ggo.in)

ADD_EXECUTABLE(gPMC gPMC.cpp ${gPMC_GGO_C})

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/include")

MESSAGE(STATUS "CMake of gPMC was found :: ${PROJECT_SOURCE_DIR}")

TARGET_LINK_LIBRARIES(gPMC
  ${GPMC_LIBRARIES}
  ${ITK_LIBRARIES}
  OpenCL::OpenCL
)

INSTALL( TARGETS  gPMC DESTINATION bin )
SET(GPMC_DLLS
  "${GPMC_LIB_DIR}/goPMC.dll"
  "${GPMC_LIB_DIR}/dcmtk.dll"
  "${GPMC_LIB_DIR}/libDicomRT.dll"
)

INSTALL( FILES ${GPMC_DLLS} DESTINATION bin )

FILE( GLOB LUT_FILES "${PROJECT_SOURCE_DIR}/input/*.bin" )
INSTALL( FILES ${LUT_FILES} DESTINATION lut )
