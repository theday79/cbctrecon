CMAKE_MINIMUM_REQUIRED (VERSION 3.8) # Yes, less won't do.

PROJECT(CbctReconApp)
# set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

IF(APPLE)
  SET(CMAKE_INSTALL_RPATH ${CMAKE_CURRENT_BINARY_DIR}/../bin)
ENDIF(APPLE)

FIND_PACKAGE(DCMTK)
FIND_PACKAGE(ITK)
FIND_PACKAGE(RTK)
FIND_PACKAGE(Plastimatch)
# FIND_PACKAGE(TBB)

# Find BZIP2 (required for Varian HNC)
OPTION(USE_BZIP2 "RTK has enabled hnc" OFF) # Hnc not implemented yet
IF(USE_BZIP2)
  FIND_PACKAGE(BZip2)
  INCLUDE_DIRECTORIES(${BZIP2_INCLUDE_DIR})
ENDIF(USE_BZIP2)

OPTION(USE_CUDA "Toggle CUDA on or off" ON)
OPTION(USE_OPENCL "Toggle OpenCL on or off" OFF) # OFF because currently broken
OPTION(USE_OPENMP "Toggle OpenMP on or off" ON)
OPTION(USE_GPMC "Toggle gPMC on or off" ON)

include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
  set(N_THREADS ${N} CACHE STRING "Number of available threads including hyperthreads")
else()
  set(N_THREADS 4 CACHE STRING "Number of available threads including hyperthreads")
endif()
ADD_DEFINITIONS(-DN_THREADS=${N_THREADS})
add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON)

IF(USE_GPMC)
  message (STATUS "Compiling with GPMC option")
  ADD_DEFINITIONS(-DUSE_GPMC=TRUE)
ENDIF(USE_GPMC)
IF(USE_CUDA)
  message (STATUS "Compiling with CUDA option")
  ADD_DEFINITIONS(-DUSE_CUDA=TRUE)
ENDIF(USE_CUDA)
IF(USE_OPENCL)
  message (STATUS "Compiling with OPENCL option")
  ADD_DEFINITIONS(-DUSE_OPENCL=TRUE)
ENDIF(USE_OPENCL)
IF(USE_CLFFT)
  message (STATUS "Compiling with clFFT option")
  ADD_DEFINITIONS(-DUSE_CLFFT=TRUE)
ENDIF(USE_CLFFT)

IF (OpenCL_FOUND)
  ADD_DEFINITIONS(-DUSE_OPENCL_PLM=TRUE)
  ADD_DEFINITIONS(-DUSE_OPENCL_RTK=TRUE)
  IF(USE_CLFFT)
    FIND_PACKAGE(clFFT REQUIRED)
  ENDIF(USE_CLFFT)
ENDIF (OpenCL_FOUND)

IF(USE_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF(OpenMP_FOUND)
    message (STATUS "Compiling with OPENMP option")
    ADD_DEFINITIONS(-DUSE_OPENMP=TRUE)
  ENDIF()
ENDIF(USE_OPENMP)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

############## Qt start wrapping ##############
FIND_PACKAGE(Qt5 COMPONENTS Core Gui OpenGL PrintSupport Widgets REQUIRED)

SET(CbctRecon_SOURCES
	main.cpp
	cbctrecon.cpp
	YK16GrayImage.cxx
	AG17RGBAImage.cxx
	qyklabel.cpp
	qcustomplot.cpp
	DlgRegistration.cxx
	StructureSet.cxx
	DlgExternalCommand.cxx
    rtkOpenCLUtilities.cxx
    rtkOpenCLFDKBackProjectionImageFilter.cxx
	OpenCLFFTFilter.cxx
	CMakeSetup.rc
)
SET(CbctRecon_HEADERS
	cbctrecon.h
	YK16GrayImage.h
	AG17RGBAImage.h
	qyklabel.h
	qcustomplot.h
	DlgRegistration.h
	StructureSet.h
	DlgExternalCommand.h
	rtkOpenCLFDKConeBeamReconstructionFilter.h
	rtkOpenCLFDKConeBeamReconstructionFilter.hxx
	rtkOpenCLFDKBackProjectionImageFilter.h
	rtkOpenCLFFTRampImageFilter.h
	rtkOpenCLFFTConvolutionImageFilter.h
	rtkOpenCLFFTConvolutionImageFilter.hxx
	itkOpenCLImageToImageFilter.h
	itkOpenCLImageToImageFilter.hxx
    rtkOpenCLUtilities.h
	OpenCLFFTFilter.h
)
SET(CbctRecon_FORMS
	cbctrecon.ui
	DlgRegistration.ui
	DlgHistogram.ui
	DlgExternalCommand.ui
)
QT5_WRAP_UI(CbctRecon_FORMS_HEADERS ${CbctRecon_FORMS})

# Executable
ADD_EXECUTABLE(CbctRecon
  ${CbctRecon_SOURCES}
  ${CbctRecon_HEADERS}
  ${CbctRecon_FORMS_HEADERS}
)

SET(QT_USE_QTCORE TRUE)
SET(QT_USE_QTGUI TRUE)
SET(QT_USE_QTOPENGL TRUE)
############## Qt end wrapping ##############

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE(${ITK_USE_FILE})
INCLUDE(${RTK_USE_FILE})

IF(WIN32)
  FUNCTION(PREPEND var prefix)
    SET(listVar "")
    FOREACH(f ${ARGN})
      GET_FILENAME_COMPONENT(extension ${f} EXT)
      GET_FILENAME_COMPONENT(file_name ${f} NAME)
      MESSAGE(STATUS ${file_name})
      IF(${file_name} MATCHES "cudart.lib")
        MESSAGE(STATUS "cudart.lib was removed from linker")
      ELSEIF(extension MATCHES ".lib")
        LIST(APPEND listVar "${f}")
      ELSE()
        LIST(APPEND listVar "${prefix}/${f}.lib")
      ENDIF()
    ENDFOREACH(f)
    SET(${var} "${listVar}" PARENT_SCOPE)
  ENDFUNCTION(PREPEND)
ELSE()
  FUNCTION(PREPEND var prefix)
    SET(listVar "")
    FOREACH(f ${ARGN})
      GET_FILENAME_COMPONENT(extension ${f} EXT)
      GET_FILENAME_COMPONENT(file_name ${f} NAME)
      MESSAGE(STATUS ${file_name})
      LIST(APPEND listVar "${prefix}/lib${f}.a")
    ENDFOREACH(f)
    SET(${var} "${listVar}" PARENT_SCOPE)
  ENDFUNCTION(PREPEND)
ENDIF()

PREPEND(RTK_LIBS_FULLPATH "${RTK_DIR}/bin/${CMAKE_CFG_INTDIR}" ${RTK_LIBRARIES})


set(PLASTIMATCH_INCLUDE_DIRS
    ${PLASTIMATCH_INCLUDE_DIRS}
    "${Plastimatch_SRC}/src/"
    "${Plastimatch_SRC}/src/plastimatch/"
    "${Plastimatch_SRC}/src/plastimatch/util"
    "${Plastimatch_SRC}/src/plastimatch/dose"
    "${Plastimatch_SRC}/src/plastimatch/sys"
    "${Plastimatch_SRC}/src/plastimatch/base"
    "${Plastimatch_SRC}/src/plastimatch/register"
    "${Plastimatch_SRC}/src/plastimatch/reconstruct"
    "${Plastimatch_SRC}/src/plastimatch/opencl"
    "${Plastimatch_SRC}/src/plastimatch/segment"
    "${Plastimatch_BUILD_DIR}/")

#this is for "C++ General-> Add below to Additional Include Directories"
INCLUDE_DIRECTORIES(
  ${ITK_INCLUDE_DIRS}
  ${RTK_INCLUDE_DIRS}
  ${PLASTIMATCH_INCLUDE_DIRS}
  ${DCMTK_INCLUDE_DIR}
)
IF(USE_CLFFT)
  INCLUDE_DIRECTORIES(${CLFFT_INCLUDE_DIRS})
ELSE()
  SET(CLFFT_LIBRARIES "")
ENDIF()

message (STATUS "ITK_LIBRARY_DIRS = ${ITK_LIBRARIES}")
message (STATUS "RTK_LIBRARY_DIRS = ${RTK_LIBRARIES}")
message (STATUS "DCMTK_LIBRARY_DIRS = ${DCMTK_LIBRARIES}")
message (STATUS "PLAS_LIBRARY_DIRS = ${PLASTIMATCH_LIBRARIES}")
message (STATUS "BZIP_LIBRARY_DIRS = ${BZIP2_LIBRARIES}")
message (STATUS "CLFFT_LIBRARY_DIRS = ${CLFFT_LIBRARIES}")

if(WIN32)
  SET(WIN_LIB netapi32.lib)
endif(WIN32)

set(QT_LIBRARIES)
list(APPEND QT_LIBRARIES
  "Qt5::Core"
  "Qt5::Gui"
  "Qt5::OpenGL"
  "Qt5::PrintSupport"
  "Qt5::Widgets"
  )

#this is for "Linker->Additional dependencies"
TARGET_LINK_LIBRARIES(CbctRecon
${WIN_LIB}
${ITK_LIBRARIES}
${RTK_LIBS_FULLPATH}
${DCMTK_LIBRARIES}
${PLASTIMATCH_LIBRARIES}
${BZIP2_LIBRARIES}
${CLFFT_LIBRARIES}
Qt5::Core
Qt5::Gui
Qt5::OpenGL
Qt5::PrintSupport
Qt5::Widgets
${TBB_LIBRARIES}
OpenCL::OpenCL
)

############## INSTALL ##############
## Using install/bin as storage for all dependencies before copy to INSTALL_PREFIX/bin ##

if($<CONFIG:Debug>)
  get_target_property(WIDGETDLL Qt5::Widgets      IMPORTED_LOCATION_DEBUG)
  get_target_property(GUIDLL    Qt5::Gui          IMPORTED_LOCATION_DEBUG)
  get_target_property(OPENGLDLL Qt5::OpenGL       IMPORTED_LOCATION_DEBUG)
  get_target_property(COREDLL   Qt5::Core         IMPORTED_LOCATION_DEBUG)
  get_target_property(COREDLL   Qt5::PrintSupport IMPORTED_LOCATION_DEBUG)
else($<CONFIG:Debug>)
  get_target_property(WIDGETDLL Qt5::Widgets      IMPORTED_LOCATION_RELEASE)
  get_target_property(GUIDLL    Qt5::Gui          IMPORTED_LOCATION_RELEASE)
  get_target_property(OPENGLDLL Qt5::OpenGL       IMPORTED_LOCATION_RELEASE)
  get_target_property(COREDLL   Qt5::Core         IMPORTED_LOCATION_RELEASE)
  get_target_property(COREDLL   Qt5::PrintSupport IMPORTED_LOCATION_RELEASE)
endif($<CONFIG:Debug>)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../bin)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../lut)

if(WIN32)
add_custom_command(
  TARGET CbctRecon POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_FILE:Qt5::Widgets>
    $<TARGET_FILE:Qt5::Gui>
    $<TARGET_FILE:Qt5::OpenGL>
    $<TARGET_FILE:Qt5::Core>
    $<TARGET_FILE:Qt5::PrintSupport>
    ${CMAKE_CURRENT_SOURCE_DIR}/fdk_opencl.cl
    ${CMAKE_CURRENT_BINARY_DIR}/../bin
)
elseif(UNIX)
  foreach(target ${QT_LIBRARIES})
    get_target_property(type ${target} TYPE)
    if(NOT type STREQUAL "SHARED_LIBRARY")
      continue()
    endif()
    get_property(location TARGET ${target} PROPERTY LOCATION_RELEASE)
    # Install .so and versioned .so.x.y
    get_filename_component(QT_LIB_DIR_tmp ${location} PATH)
    get_filename_component(QT_LIB_NAME_tmp ${location} NAME)
    string(REPLACE ".${_qt_version}" "" QT_LIB_NAME_tmp ${QT_LIB_NAME_tmp})
    message(STATUS "Attmpting installing ${QT_LIB_NAME_tmp} from dir: ${QT_LIB_DIR_tmp}")
    install(DIRECTORY ${QT_LIB_DIR_tmp}/
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin/${QT_LIB_NAME_tmp}.framework/ COMPONENT Runtime
      FILES_MATCHING PATTERN "${QT_LIB_NAME_tmp}*"
      PATTERN "${QT_LIB_NAME_tmp}*.debug" EXCLUDE)
  endforeach()
endif()

IF(USE_CLFFT AND WIN32)
add_custom_command(
  TARGET CbctRecon POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CLFFT_ROOT_DIR}/bin/clFFT.dll
    ${CMAKE_CURRENT_BINARY_DIR}/../bin
)
ENDIF()


INSTALL( TARGETS CbctRecon DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin )
