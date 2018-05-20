CMAKE_MINIMUM_REQUIRED (VERSION 3.8) # Yes, less won't do.

PROJECT(App)
# set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

IF(APPLE)
  SET(CMAKE_INSTALL_RPATH bin) # ${CMAKE_CURRENT_BINARY_DIR}/../bin
ENDIF(APPLE)

FIND_PACKAGE(DCMTK)
FIND_PACKAGE(ITK REQUIRED)
FIND_PACKAGE(Plastimatch REQUIRED)
if(NOT RTK_AS_ITK_EXTERNAL)
  FIND_PACKAGE(RTK REQUIRED)
endif()

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
	PlmWrapper.cxx
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
	PlmWrapper.h
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

link_directories(${ITK_LIBRARY_DIRS})

# Executable
ADD_EXECUTABLE(CbctRecon
  ${CbctRecon_SOURCES}
  ${CbctRecon_HEADERS}
  ${CbctRecon_FORMS_HEADERS}
)

# set_target_properties(CbctRecon PROPERTIES POSITION_INDEPENDENT_CODE ON) #lib only property?

# set_property(TARGET CbctRecon
#              PROPERTY CUDA_SEPARABLE_COMPILATION ${CUDA_SEPARABLE_COMPILATION})

SET(QT_USE_QTCORE TRUE)
SET(QT_USE_QTGUI TRUE)
SET(QT_USE_QTOPENGL TRUE)
############## Qt end wrapping ##############

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# INCLUDE(${ITK_USE_FILE})
# INCLUDE(${RTK_USE_FILE})

IF(RTK_AS_ITK_EXTERNAL)
  set(RTK_LIBS_FULLPATH ${RTK_LIBRARIES})
ELSE()
  IF(WIN32)
    FUNCTION(PREPEND var prefix)
      SET(listVar "")
      FOREACH(f ${ARGN})
        GET_FILENAME_COMPONENT(extension ${f} EXT)
        GET_FILENAME_COMPONENT(file_name ${f} NAME)
        # MESSAGE(STATUS ${file_name})
        #IF(${file_name} MATCHES "cudart.lib")
        #  MESSAGE(STATUS "cudart.lib was removed from linker")
        #ELSE
        IF(extension MATCHES ".lib")
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
ENDIF()

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
message (STATUS "RTK_LIBRARIES = ${RTK_LIBRARIES}")
message (STATUS "RTK_LIBRARY_DIRS = ${RTK_LIBRARY_DIRS}")
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
  set(QT_DLLS
    $<TARGET_FILE:Qt5::Widgets>
    $<TARGET_FILE:Qt5::Gui>
    $<TARGET_FILE:Qt5::OpenGL>
    $<TARGET_FILE:Qt5::Core>
    $<TARGET_FILE:Qt5::PrintSupport>
  )
  install(FILES ${QT_DLLS} DESTINATION bin)
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



INSTALL( TARGETS CbctRecon DESTINATION bin )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/fdk_opencl.cl DESTINATION bin)
# ${CMAKE_CURRENT_BINARY_DIR}/../bin

##-----------------------------------------------------------------------------
##  Additional install files
##-----------------------------------------------------------------------------
if (WIN32)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/WindeployQt.cmake)
  if($<CONFIG:Debug>)
    windeployqt(CbctRecon "debug" "${CMAKE_CURRENT_BINARY_DIR}/QtLibs")
  else()
    windeployqt(CbctRecon "release" "${CMAKE_CURRENT_BINARY_DIR}/QtLibs")
  endif()

  file (GLOB DLLFILES "${CMAKE_CURRENT_BINARY_DIR}/QtLibs/*.dll")
  install (FILES ${DLLFILES} DESTINATION bin)
  set (QT_SUBDIRS "iconengines" "imageformats" "platforms" "styles" "translations")
  foreach(subdir ${QT_SUBDIRS})
    install(DIRECTORY "$<TARGET_FILE_DIR:CbctRecon>/${subdir}" DESTINATION bin)
  endforeach()

  if (OPENMP_FOUND)
    set(CMAKE_INSTALL_OPENMP_LIBRARIES ON)
    # install(FILES "${MSVC_REDIST_DIR}/x64/Microsoft.VC140.OpenMP/vcomp140.dll" DESTINATION bin COMPONENT Libraries)
    # install(FILES "C:/Program Files (x86)/Microsoft Visual Studio 9.0/VC/redist/x86/Microsoft.VC90.OPENMP/Microsoft.VC90.OpenMP.manifest" DESTINATION bin COMPONENT Libraries)
  endif(OPENMP_FOUND)
  include(InstallRequiredSystemLibraries)

	if (CUDA_FOUND)
	  set (DLL_PATTERNS "cufft64*.dll" "cublas64*.dll" "cudart64*.dll")
	  foreach (p ${DLL_PATTERNS})
		file (GLOB DLLFILES "${CUDA_TOOLKIT_ROOT_DIR}/bin/${p}")
		if (DLLFILES)
		  install (FILES ${DLLFILES} DESTINATION bin)
		endif ()
	  endforeach ()
	endif ()

	if (ITK_FOUND)
    # ITK*.dll might not include RTK_AS_ITK_EXTERNAL libs
    set (DLL_PATTERNS "*.dll")
	  foreach (p ${DLL_PATTERNS})
      file (GLOB DLLFILES "${ITK_DIR}/bin/release/${p}")
      if (DLLFILES)
        install (FILES ${DLLFILES} DESTINATION bin)
      endif ()
	  endforeach ()
	endif ()

	if (RTK_FOUND)
	  install (FILES ${RTK_DIR}/bin/release/RTK.dll DESTINATION bin)
	  install (FILES ${RTK_DIR}/bin/release/rtkcuda.dll DESTINATION bin)
	  install (FILES ${RTK_DIR}/bin/release/ITKCudaCommon.dll DESTINATION bin)
	  install (FILES ${RTK_DIR}/bin/release/lpsolve55.dll DESTINATION bin) #is 55 a version number?
	endif ()

	if (Plastimatch_FOUND)
	  set (DLL_PATTERNS "plm*.dll" "")
	  foreach (p ${DLL_PATTERNS})
      file (GLOB DLLFILES "${Plastimatch_BUILD_DIR}/release/${p}")
      if (DLLFILES)
        install (FILES ${DLLFILES} DESTINATION bin)
      endif ()
	  endforeach ()
	endif ()

	IF(USE_CLFFT)
	  install (FILES ${CLFFT_ROOT_DIR}/bin/clFFT.dll DESTINATION bin)
	ENDIF()

	#Add FFT dlls to Install
	if (FFTW_FOUND)
	  if (EXISTS "${FFTW_DIR}/libfftw3-3.dll")
		  install (FILES "${FFTW_DIR}/libfftw3-3.dll" DESTINATION bin)
	  endif ()
	  if (EXISTS "${FFTW_DIR}/libfftw3f-3.dll")
		  install (FILES "${FFTW_DIR}/libfftw3f-3.dll" DESTINATION bin)
	  endif ()
	endif ()
endif()

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
