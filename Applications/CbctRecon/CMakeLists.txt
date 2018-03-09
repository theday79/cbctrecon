CMAKE_MINIMUM_REQUIRED (VERSION 3.8) # Yes, less won't do.

PROJECT(App)
# set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

IF(APPLE)
  SET(CMAKE_INSTALL_RPATH bin) # ${CMAKE_CURRENT_BINARY_DIR}/../bin
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
OPTION(USE_LOWPASS_FFT "Gauss FFT as low-pass filter for scatter smoothing" ON)

# Packaging
option (CbctRecon_PACKAGE_32BIT
    "Set this when building 32-bit packages on a 64-bit machine" OFF)
option (CbctRecon_PACKAGE_NSIS "Set to ON when packaging binaries with NSIS" OFF)
option (CbctRecon_PACKAGE_WIX "Set to ON when packaging binaries with WIX" ON)

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

IF(USE_LOWPASS_FFT)
  message (STATUS "Compiling with Gauss FFT as low-pass filter for scatter smoothing option")
  ADD_DEFINITIONS(-DLOWPASS_FFT=TRUE)
ENDIF(USE_LOWPASS_FFT)

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
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF(USE_OPENMP)


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
	DlgHistogram.ui # Kept for example, no longer used
	DlgExternalCommand.ui
)
QT5_WRAP_UI(CbctRecon_FORMS_HEADERS ${CbctRecon_FORMS})

set_property(SOURCE ui_cbctrecon.h PROPERTY SKIP_AUTOMOC ON)
set_property(SOURCE ui_DlgRegistration.h PROPERTY SKIP_AUTOMOC ON)
set_property(SOURCE ui_DlgHistogram.h PROPERTY SKIP_AUTOMOC ON)
set_property(SOURCE ui_DlgExternalCommand.h PROPERTY SKIP_AUTOMOC ON)

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

file(MAKE_DIRECTORY bin) # ${CMAKE_CURRENT_BINARY_DIR}/../bin
file(MAKE_DIRECTORY lut) # ${CMAKE_CURRENT_BINARY_DIR}/../lut

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
      DESTINATION bin/${QT_LIB_NAME_tmp}.framework/ COMPONENT Runtime # ${CMAKE_CURRENT_BINARY_DIR}/../
      FILES_MATCHING PATTERN "${QT_LIB_NAME_tmp}*"
      PATTERN "${QT_LIB_NAME_tmp}*.debug" EXCLUDE)
  endforeach()
endif()

IF(USE_CLFFT AND WIN32)
add_custom_command(
  TARGET CbctRecon POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CLFFT_ROOT_DIR}/bin/clFFT.dll
    bin # ${CMAKE_CURRENT_BINARY_DIR}/../bin
)
ENDIF()


INSTALL( TARGETS CbctRecon DESTINATION bin )
# ${CMAKE_CURRENT_BINARY_DIR}/../bin


##-----------------------------------------------------------------------------
##  Packaging - Copied with modification from Plastimatch
##-----------------------------------------------------------------------------
##  Simple version: "make package_source" to create tarball
##
##  Advanced version: to generate...
##   * zip    : cpack -G ZIP <source tree>
##   * tar.bz2: cpack -G TBZ2 <source tree>
##
##    ...from within your build directory
##-----------------------------------------------------------------------------
# Choose generator
if (CbctRecon_PACKAGE_WIX)
  set (CPACK_GENERATOR "WIX")
elseif (CbctRecon_PACKAGE_NSIS)
  set (CPACK_GENERATOR "NSIS")
else ()
  set (CPACK_GENERATOR "ZIP")
endif ()
set (CPACK_SOURCE_GENERATOR "TBZ2")

# For WiX install, and possibly for debian -dev targets,
# we should properly set up the install components.
# Until then, let's just kill it
set (CPACK_MONOLITHIC_INSTALL 1)

# General stuff
set (CPACK_INSTALL_PREFIX "")
set (CPACK_PACKAGE_CONTACT "cbctrecon@googlegroups.com")
set (CPACK_PACKAGE_NAME "cbctrecon")
set (CPACK_PACKAGE_INSTALL_DIRECTORY "CbctRecon")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "CbctRecon - A Medical Imaging Application")
set (CPACK_PACKAGE_VENDOR "CbctRecon Dev Team")
set (CPACK_PACKAGE_VERSION_MAJOR "${CbctRecon_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${CbctRecon_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${CbctRecon_VERSION_PATCH}")
set (CPACK_PACKAGE_VERSION_TWEAK "${CbctRecon_VERSION_TWEAK}")
set (CPACK_SOURCE_IGNORE_FILES
  "/\\\\.git"
  "/\\\\.svn"
  "~$"
)

# NSIS stuff
set(CPACK_NSIS_INSTALLED_ICON_NAME "Icon1.ico")
set (CPACK_NSIS_HELP_LINK "https://github.com/agravgaard/cbctrecon")
set (CPACK_NSIS_URL_INFO_ABOUT "https://github.com/agravgaard/cbctrecon")
set (CPACK_NSIS_CONTACT "andreasga22@gmail.com")
set (CbctRecon_NSIS_VERSION_STRING "${CbctRecon_VERSION_MAJOR}.${CbctRecon_VERSION_MINOR}.${CbctRecon_VERSION_PATCH}")

# WIX stuff
set (CbctRecon_WINDOWS_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extra/windows-install")
set (CPACK_WIX_HELP_LINK "https://github.com/agravgaard/cbctrecon")
set (CPACK_WIX_UPGRADE_GUID "43DF18E5-D982-413F-834C-2865B977C145")
set (CPACK_WIX_LICENSE_RTF  "${CbctRecon_WINDOWS_INSTALL_DIR}/License.rtf")
set (CPACK_WIX_UI_BANNER  "${CbctRecon_WINDOWS_INSTALL_DIR}/cbctrecon_logo.jpg")
set (CPACK_WIX_UI_DIALOG  "${CbctRecon_WINDOWS_INSTALL_DIR}/cbctrecon_logo_GUI.jpg")

set (CPACK_WIX_EXTRA_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/extra/wix/cmake_extra_dialog.wxs"
    "${CMAKE_CURRENT_SOURCE_DIR}/extra/wix/install_dir.wxs"
    "${CMAKE_CURRENT_SOURCE_DIR}/extra/wix/start_menu.wxs"
    )
set (CPACK_WIX_PATCH_FILE
    "${CMAKE_CURRENT_SOURCE_DIR}/extra/wix/patch_path_env.xml"
#    "${CMAKE_CURRENT_SOURCE_DIR}/extra/wix/patch_start_menu.xml"
    )
set (CPACK_WIX_UI_REF "CMakeUI_InstallDir")
if (CPACK_GENERATOR STREQUAL "WIX")
  set (CPACK_PACKAGE_NAME "CbctRecon")
  if (CPACK_PACKAGE_VERSION_TWEAK AND CPACK_PACKAGE_VERSION_TWEAK GREATER "0")
    set (CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.${CPACK_PACKAGE_VERSION_TWEAK}")
    if (CMAKE_SIZEOF_VOID_P EQUAL 8)
      set (CPACK_SYSTEM_NAME "win64")
    else ()
      set (CPACK_SYSTEM_NAME "win32")
    endif ()
    set (CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME}")
  endif ()
  set_property (INSTALL "bin/CbctRecon.exe"
    PROPERTY CPACK_START_MENU_SHORTCUTS "CbctRecon Application")
endif ()

# ZIP stuff
if (CPACK_GENERATOR STREQUAL "ZIP")
  set (CPACK_SET_DESTDIR ON)
endif ()

include (CPack)
