CMAKE_MINIMUM_REQUIRED (VERSION 3.8)


# This project is designed to be built outside the RTK source tree.


#test
PROJECT(CbctRecon)
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
if(POLICY CMP0020)
  cmake_policy(SET CMP0020 OLD)
endif()
# Find the RTK libraries and includes


SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

FIND_PACKAGE(RTK REQUIRED)
INCLUDE(${RTK_USE_FILE})

FIND_PACKAGE(CUDA)

FIND_PACKAGE(OpenCL REQUIRED)

get_filename_component(MY_PARENT_DIR ${CUDA_LIBRARIES} PATH)
set (CUDA_FFT_LIB ${MY_PARENT_DIR}/cufft.lib)

find_package(Qt5 COMPONENTS Core Gui OpenGL PrintSupport Widgets REQUIRED)


find_package (Plastimatch)
set(PLASTIMATCH_INCLUDE_DIRS
    ${PLASTIMATCH_INCLUDE_DIRS}
    "${Plastimatch_SRC}/src/"
    "${Plastimatch_SRC}/src/plastimatch/"
    "${Plastimatch_SRC}/src/plastimatch/util"
    "${Plastimatch_SRC}/src/plastimatch/dose"
    "${Plastimatch_SRC}/src/plastimatch/sys"
    "${Plastimatch_SRC}/src/plastimatch/base"
    "${Plastimatch_SRC}/src/plastimatch/register"
    "${Plastimatch_SRC}/src/plastimatch/segment"
    "${Plastimatch_BUILD_DIR}/")

message (STATUS "Plastimatch_FOUND = ${Plastimatch_FOUND}")
message (STATUS "Plastimatch_DIR = ${Plastimatch_DIR}")
message (STATUS "Plastimatch_SRC = ${Plastimatch_SRC}")
message (STATUS "Plastimatch_BUILD_DIR = ${Plastimatch_BUILD_DIR}")

message (STATUS "CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}")
message (STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}")
message (STATUS "MY_PARENT_DIR = ${MY_PARENT_DIR}")

message (STATUS "PLASTIMATCH_LIBRARIES = ${PLASTIMATCH_LIBRARIES}")
message (STATUS "PLASTIMATCH_INCLUDE_DIRS = ${PLASTIMATCH_INCLUDE_DIRS}")

#Find ITk folder
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
   INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

#Find DCMTK

find_package(DCMTK REQUIRED)
LINK_DIRECTORIES("${DCMTK_DIR}/../lib")

#works But A lot of lInk error (detection...

MESSAGE(STATUS "DCMTK_FOUND = ${DCMTK_FOUND}")
MESSAGE(STATUS "DCMTK_DIR = ${DCMTK_DIR}")

## JUST TEMPORARY WHILE USING THE OLD OPENCL VERSION OF RTK ##
# Find BZIP2 (required for Varian HNC)
OPTION(USE_BZIP2 "RTK has enabled hnc" OFF)
IF(USE_BZIP2)
  FIND_PACKAGE(BZip2)
  IF(BZIP2_FOUND)
    INCLUDE_DIRECTORIES(${BZIP2_INCLUDE_DIR})
  ELSE(BZIP2_FOUND)
    MESSAGE(FATAL_ERROR "Cannot build without bzip2.")
  ENDIF(BZIP2_FOUND)
ENDIF(USE_BZIP2)

OPTION(USE_GPMC "Enable gpu Proton Monte Carlo (gPMC) dose re-calculation" ON)
IF(USE_GPMC)
  include(ExternalProject)
  ExternalProject_Add(gPMC ## Yes, external is needed because otherwise dlls will interfere with local DCMTK
	DOWNLOAD_COMMAND ""    ## Also, DCMTK must be installed in a path not called %PROGRAMFILES%/DCMTK !!!
	SOURCE_DIR "${PROJECT_SOURCE_DIR}/gPMC"
	# CONFIGURE_COMMAND ""
	CMAKE_GENERATOR "Visual Studio 12 2013 Win64"
	CMAKE_GENERATOR_TOOLSET "v120" ## Visual Studio 2013
	CONFIGURATIONS "Debug" ## Release has not worked yet, probably some "if debug" directives in goPMC or dependencies
	INSTALL_DIR "${CMAKE_BINARY_DIR}/install"
	CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/install"
  )
  ADD_DEFINITIONS(-DUSE_GPMC=TRUE)
ENDIF(USE_GPMC)

SET(CbctRecon_SOURCES
	main.cpp
	cbctrecon.cpp
	YK16GrayImage.cxx
	qyklabel.cpp
	qcustomplot.cpp
	DlgRegistration.cxx
	DlgExternalCommand.cxx
)

SET(CbctRecon_HEADERS
	cbctrecon.h
	YK16GrayImage.h
	qyklabel.h
	qcustomplot.h
	DlgRegistration.h
	DlgExternalCommand.h
)

SET(CbctRecon_FORMS
	cbctrecon.ui
	DlgRegistration.ui
	DlgExternalCommand.ui
)



QT5_WRAP_UI(CbctRecon_FORMS_HEADERS ${CbctRecon_FORMS})


# Executable
ADD_EXECUTABLE(CbctRecon
  ${CbctRecon_SOURCES}
  ${CbctRecon_HEADERS}
  ${CbctRecon_FORMS_HEADERS}
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})


message (STATUS "PLASTIMATCH_SOURCE_DIRECTORY = ${PLASTIMATCH_INCLUDE_DIRS}")
message (STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")


#this is for "C++ General-> Add below to Additional Include Directories"
INCLUDE_DIRECTORIES(${PLASTIMATCH_INCLUDE_DIRS}
  ${DCMTK_INCLUDE_DIR}
)

SET(QT_USE_QTCORE TRUE)
SET(QT_USE_QTGUI TRUE)
SET(QT_USE_QTOPENGL TRUE)

message (STATUS "RTK_USE_FILE = ${RTK_USE_FILE}")
message (STATUS "RTK_LIBRARY_DIRS = ${RTK_LIBRARY_DIRS}")

#this is for "Linker->Additional dependencies"
TARGET_LINK_LIBRARIES(CbctRecon
netapi32.lib
${ITK_LIBRARIES}
${RTK_LIBRARIES}
${GPMC_LIBRARIES}
${PLASTIMATCH_LIBRARIES}
${BZIP2_LIBRARIES}
OpenCL::OpenCL
Qt5::Core
Qt5::Gui
Qt5::OpenGL
Qt5::PrintSupport
Qt5::Widgets
)

## INSTALL ## Using install/bin as storage for all dependencies before copy to INSTALL_PREFIX/bin ##

if($<CONFIG:Debug>)
  get_target_property(WIDGETDLL Qt5::Widgets IMPORTED_LOCATION_DEBUG)
else($<CONFIG:Debug>)
  get_target_property(WIDGETDLL Qt5::Widgets      IMPORTED_LOCATION_RELEASE)
  get_target_property(GUIDLL    Qt5::Gui          IMPORTED_LOCATION_RELEASE)
  get_target_property(COREDLL   Qt5::Core         IMPORTED_LOCATION_RELEASE)
  get_target_property(COREDLL   Qt5::PrintSupport IMPORTED_LOCATION_RELEASE)
endif($<CONFIG:Debug>)

add_custom_command(
  TARGET CbctRecon POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    $<TARGET_FILE:Qt5::Widgets>
    $<TARGET_FILE:Qt5::Gui>
    $<TARGET_FILE:Qt5::Core>
    $<TARGET_FILE:Qt5::PrintSupport>
    install/bin #$<TARGET_FILE_DIR:CbctRecon>
)

INSTALL( TARGETS CbctRecon DESTINATION bin )

FILE( GLOB DLLFILES "${CMAKE_CURRENT_BINARY_DIR}/install/bin/*" ) # also gPMC.exe
INSTALL( FILES ${DLLFILES} DESTINATION bin )
