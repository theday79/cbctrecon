CMAKE_MINIMUM_REQUIRED (VERSION 3.10) # Yes, less won't do.

set(CMAKE_INCLUDE_CURRENT_DIR ON)

#-----------------------------------------------------------------------------
# BASIC CMAKE
set(CbctRecon_VERSION_MAJOR "1") # bump at major change -> does NOT mean release version
set(CbctRecon_VERSION_MINOR "4") # bump at functionality change
set(CbctRecon_VERSION_PATCH "0") # bump at minor change

project(CbctRecon
  VERSION "${CbctRecon_VERSION_MAJOR}.${CbctRecon_VERSION_MINOR}.${CbctRecon_VERSION_PATCH}"
  LANGUAGES CXX
  )

SET(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ standard specification")

set(CbctRecon_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
#-----------------------------------------------------------------------------


OPTION(USE_CUDA "Toggle CUDA on or off" OFF) # OFF by default for multiple reasons
OPTION(USE_OPENMP "Toggle OpenMP on or off" ON)
OPTION(USE_CLFFT "Toggle CLFFT (in development) on or off" OFF) # OFF by default because it doesn't yet work
OPTION(USE_GPMC "Toggle gPMC on or off" OFF) # OFF by default because you can't expect developers to have access to gPMC binaries AND use VS2013 Debug when Release is much faster for everything else.
OPTION(USE_LOWPASS_FFT "Gauss FFT as low-pass filter for scatter smoothing" ON)
OPTION(CBCTRECON_CLANG_COMPLETE "Generate .clang_complete for vim" ON)
OPTION(CBCTRECON_BUILD_TESTS "Build tests for CbctRecon" ON)

OPTION(BUILD_SHARED_LIBS "Build shared libs in ITK, Plastimatch and RTK (for redistribution)" ON)
if(NOT WIN32)
  set(BUILD_SHARED_LIBS OFF)
  OPTION(USE_TINYREFL "Enable tinyrefl static reflection tool for CbctRecon tests" ON)
else()
  OPTION(USE_TINYREFL "Enable tinyrefl static reflection tool for CbctRecon tests" OFF)
endif()
OPTION(CBCTRECON_BUILD_SHARED_LIBS "Build shared libraries for CbctRecon" ${BUILD_SHARED_LIBS})

OPTION(BUILD_TESTING "Build testing in ITK, RTK and Plastimatch" OFF)
# BUILD_TESTING is OFF by default, because ITK and RTK has decent CI tests

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

IF(USE_GPMC)
  message (STATUS "Compiling with GPMC option")
  add_definitions(-DUSE_GPMC=TRUE)
ENDIF(USE_GPMC)

IF(USE_LOWPASS_FFT)
  message (STATUS "Compiling with Gauss FFT as low-pass filter for scatter smoothing option")
  ADD_DEFINITIONS(-DLOWPASS_FFT=TRUE)
ENDIF(USE_LOWPASS_FFT)

IF(USE_TINYREFL)
  message (STATUS "Compiling with tinyrefl enabled")
  ADD_DEFINITIONS(-DUSE_TINYREFL=TRUE)
ENDIF(USE_TINYREFL)

FIND_PACKAGE(OpenCL REQUIRED)
if(OpenCL_FOUND)
  add_definitions(-DUSE_OPENCL_PLM=TRUE)
  add_definitions(-DUSE_OPENCL_RTK=TRUE)
endif(OpenCL_FOUND)

IF(USE_CUDA)
  message (STATUS "Compiling with CUDA option")
  FIND_PACKAGE(CUDA_wrap)
  GET_FILENAME_COMPONENT(CUDA_TOOLKIT_ROOT_DIR "${CUDA_TOOLKIT_INCLUDE}/.." ABSOLUTE)
  FIND_PACKAGE(CUDA REQUIRED)
  add_definitions(-DUSE_CUDA=TRUE)
ENDIF(USE_CUDA)

FIND_PACKAGE(FFTW)

# set(TBB_ROOT "" CACHE PATH "Root dir of intel TBB, (one level up from lib)")
# FIND_PACKAGE(TBB) ## The world is not ready for this yet

IF(USE_OPENMP)
  FIND_PACKAGE(OpenMP REQUIRED)
  message(STATUS "Compiling with OPENMP option")
  add_definitions(-DUSE_OPENMP=TRUE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF(USE_OPENMP)

IF(USE_CLFFT)
  FIND_PACKAGE(clFFT REQUIRED)
  message(STATUS "Compiling with clFFT option")
  add_definitions(-DUSE_CLFFT=TRUE)
ENDIF(USE_CLFFT)

FIND_PACKAGE(Qt5 COMPONENTS Core Gui OpenGL PrintSupport Widgets Test REQUIRED)


#-----------------------------------------------------------------------------
# Set a default build type if none was specified
#-----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND UNIX)
  mark_as_advanced(CMAKE_BUILD_TYPE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug"
    "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
set(CMAKE_INSTALL_BINDIR "${CMAKE_BINARY_DIR}/bin" CACHE PATH "" FORCE)
set(CMAKE_INSTALL_DATAROOTDIR "${CMAKE_BINARY_DIR}/share" CACHE PATH "" FORCE)
set(CMAKE_INSTALL_INCLUDEDIR "${CMAKE_BINARY_DIR}/include" CACHE PATH "" FORCE)
set(CMAKE_INSTALL_LIBDIR "${CMAKE_BINARY_DIR}/lib64" CACHE PATH "" FORCE)
set(LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib" CACHE PATH "" FORCE)
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib" CACHE PATH "" FORCE)
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin" CACHE PATH "" FORCE)
set(EXECUTEABLE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" CACHE PATH "" FORCE)

add_subdirectory(External)

if(NOT DCMTK_FOUND)
  if(NOT USE_SYSTEM_DCMTK)
    set(DCMTK_DIR "${CMAKE_BINARY_DIR}")
  endif()
  find_package(DCMTK_wrap)
endif()

if(NOT ITK_FOUND)
  if(NOT USE_SYSTEM_ITK)
    set(ITK_DIR "${CMAKE_BINARY_DIR}/ITK-build")
  endif()
  find_package(ITK REQUIRED)
  include(${ITK_USE_FILE})
endif()

if(USE_SYSTEM_RTK)
  find_package(RTK REQUIRED)
  include(${RTK_USE_FILE})
endif()

set(PLM_TARGETS
  plmbase
  plmutil
  plmdose
  plmbase
  plmregister
  plmreconstruct
  plmopencl
  plmsegment
)

if(WIN32)
  set(WIN_LIB netapi32.lib)
endif(WIN32)

set(QT_LIBRARIES)
list(APPEND QT_LIBRARIES
  "Qt5::Core"
  "Qt5::Gui"
  "Qt5::OpenGL"
  "Qt5::PrintSupport"
  "Qt5::Widgets"
  "Qt5::Test"
  )

set(CBCTRECON_DEP_LIBS
  ${WIN_LIB}
  ${DCMTK_LIBRARIES}
  ${RTK_LIBRARIES}
  ${CLFFT_LIBRARIES}
  ${QT_LIBRARIES}
  ${TBB_LIBRARIES}
  OpenCL::OpenCL
)
add_subdirectory(Applications)

if(CBCTRECON_BUILD_TESTS)
  include(${CbctRecon_SOURCE_DIR}/cmake/ctest_setup.cmake)
  message(STATUS "Configuring CbctRecon tests...")
  add_subdirectory(Testing)
endif()

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(STATUS "Generate vim YCM config file for CbctRecon at ${CbctRecon_SOURCE_DIR} ...")
    configure_file(
        ${CbctRecon_SOURCE_DIR}/.ycm_extra_conf.py.in
        ${CbctRecon_SOURCE_DIR}/.ycm_extra_conf.py
    )
endif()
