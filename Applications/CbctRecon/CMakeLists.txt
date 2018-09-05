CMAKE_MINIMUM_REQUIRED (VERSION 3.8) # Yes, less won't do.

project(App)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC_SEARCH_PATHS
  ${CMAKE_AUTOUIC_SEARCH_PATHS}
  ${CMAKE_CURRENT_SOURCE_DIR}/ui
  )

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

if(APPLE)
  SET(CMAKE_INSTALL_RPATH bin)
endif(APPLE)

# Packaging
option(CbctRecon_PACKAGE_32BIT "Set this when building 32-bit packages on a 64-bit machine" OFF)
option(CbctRecon_PACKAGE_NSIS "Set to ON when packaging binaries with NSIS" OFF)
option(CbctRecon_PACKAGE_WIX "Set to ON when packaging binaries with WIX" ON)

add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON)

if(USE_GPMC)
  message(STATUS "Compiling with GPMC option")
  add_definitions(-DUSE_GPMC=TRUE)
endif(USE_GPMC)

if(USE_CUDA)
  message(STATUS "Compiling with CUDA option")
  add_definitions(-DUSE_CUDA=TRUE)
endif(USE_CUDA)

if(USE_LOWPASS_FFT)
  message(STATUS "Compiling with Gauss FFT as low-pass filter for scatter smoothing option")
  add_definitions(-DLOWPASS_FFT=TRUE)
endif(USE_LOWPASS_FFT)

if(OpenCL_FOUND)
  add_definitions(-DUSE_OPENCL_PLM=TRUE)
  add_definitions(-DUSE_OPENCL_RTK=TRUE)
  if(USE_CLFFT)
    FIND_PACKAGE(clFFT REQUIRED)
    message(STATUS "Compiling with clFFT option")
    add_definitions(-DUSE_CLFFT=TRUE)
  endif(USE_CLFFT)
endif(OpenCL_FOUND)

if(USE_OPENMP)
  FIND_PACKAGE(OpenMP)
  if(OpenMP_FOUND)
    message(STATUS "Compiling with OPENMP option")
    add_definitions(-DUSE_OPENMP=TRUE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif(USE_OPENMP)


############## Qt start wrapping ##############
FIND_PACKAGE(Qt5 COMPONENTS Core Gui OpenGL PrintSupport Widgets REQUIRED)

SET(CbctRecon_SOURCES
  cbctrecon.cpp
  qt-utils/YK16GrayImage.cxx
  qt-utils/AG17RGBAImage.cxx
  qt-utils/qyklabel.cpp
  qt-utils/qcustomplot.cpp
  DlgRegistration.cxx
  StructureSet.cxx
  PlmWrapper.cxx
  DlgExternalCommand.cxx
  rtkExtension/rtkOpenCLUtilities.cxx
  rtkExtension/rtkOpenCLFDKBackProjectionImageFilter.cxx
  OpenCLFFTFilter.cxx
  CMakeSetup.rc
)
SET(CbctRecon_HEADERS
  cbctrecon.h
  qt-utils/YK16GrayImage.h
  qt-utils/AG17RGBAImage.h
  qt-utils/qyklabel.h
  qt-utils/qcustomplot.h
  DlgRegistration.h
  StructureSet.h
  PlmWrapper.h
  DlgExternalCommand.h
  rtkExtension/rtkOpenCLFDKConeBeamReconstructionFilter.h
  rtkExtension/rtkOpenCLFDKConeBeamReconstructionFilter.hxx
  rtkExtension/rtkOpenCLFDKBackProjectionImageFilter.h
  rtkExtension/rtkOpenCLFFTRampImageFilter.h
  rtkExtension/rtkOpenCLFFTConvolutionImageFilter.h
  rtkExtension/rtkOpenCLFFTConvolutionImageFilter.hxx
  rtkExtension/itkOpenCLImageToImageFilter.h
  rtkExtension/itkOpenCLImageToImageFilter.hxx
  rtkExtension/rtkOpenCLUtilities.h
  OpenCLFFTFilter.h
  io.hxx
  compute.hxx
  fdk.hxx
)
SET(CbctRecon_FORMS
  ui/cbctrecon.ui
  ui/DlgRegistration.ui
  ui/DlgHistogram.ui # Kept for example, no longer used
  ui/DlgExternalCommand.ui
)

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

set(LOCAL_INCL_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}/qt-utils
  ${CMAKE_CURRENT_SOURCE_DIR}/rtkExtension
  )

if(MSVC AND NOT(MSVC_VERSION LESS 1400))
  set_source_files_properties(cbctrecon.cpp
    PROPERTIES COMPILE_FLAGS "/bigobj")
endif()

# Executable
add_executable(CbctRecon
  main.cpp
  ${CbctRecon_FORMS}
)

add_library(CbctReconLib
  ${CbctRecon_SOURCES}
  ${CbctRecon_HEADERS}
  ${CbctRecon_FORMS}
  )

set_target_properties(CbctRecon
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

SET(QT_USE_QTCORE TRUE)
SET(QT_USE_QTGUI TRUE)
SET(QT_USE_QTOPENGL TRUE)
############## Qt end wrapping ##############

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

set(Plastimatch_SRC ${CMAKE_BINARY_DIR}/Plastimatch-src/)
set(Plastimatch_BUILD_DIR ${CMAKE_BINARY_DIR}/Plastimatch-build/)

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
set(CBCTRECON_DEP_INCLUDES
  ${LOCAL_INCL_DIRS}
  ${ITK_INCLUDE_DIRS}
  ${RTK_INCLUDE_DIRS}
  ${PLASTIMATCH_INCLUDE_DIRS}
  ${DCMTK_INCLUDE_DIR}
  ${CLFFT_INCLUDE_DIRS}
)
get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(is_multi_config)
  set(GEN_UI_HEADER_DIR ${CMAKE_BINARY_DIR}/Applications/CbctRecon/CbctReconLib_autogen/include_<CONFIG>)
else()
  set(GEN_UI_HEADER_DIR ${CMAKE_BINARY_DIR}/Applications/CbctRecon/CbctReconLib_autogen/include)
endif()

target_include_directories(CbctReconLib
  PUBLIC ${CBCTRECON_DEP_INCLUDES}
  PUBLIC ${GEN_UI_HEADER_DIR}
)

target_include_directories(CbctRecon
  PUBLIC ${CBCTRECON_DEP_INCLUDES}
  PUBLIC ${GEN_UI_HEADER_DIR}
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
  )

#this is for "Linker->Additional dependencies"
set(CBCTRECON_DEP_LIBS
  ${WIN_LIB}
  ${DCMTK_LIBRARIES}
  lpsolve55 # Has to be added manually on windows?
  ${PLM_TARGETS}
  ${CLFFT_LIBRARIES}
  ${QT_LIBRARIES}
  ${TBB_LIBRARIES}
  OpenCL::OpenCL
)
target_link_libraries(CbctReconLib
  PUBLIC ${CBCTRECON_DEP_LIBS}
)
target_link_libraries(CbctRecon
  PUBLIC CbctReconLib
)
if(USE_TINYREFL)
  include(${TINYREFL_SOURCE_DIR}/tool/driver.cmake)
  tinyrefl_tool(TARGET CbctReconLib
    HEADERS
    cbctrecon.h
    COMPILE_OPTIONS
    -fPIC
    )
endif()

############## INSTALL ##############
## Using install/bin as storage for all dependencies before copy to INSTALL_PREFIX/bin ##

file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin)
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lut)

if(WIN32)
  set(QT_DLLS
    $<TARGET_FILE:Qt5::Widgets>
    $<TARGET_FILE:Qt5::Gui>
    $<TARGET_FILE:Qt5::OpenGL>
    $<TARGET_FILE:Qt5::Core>
    $<TARGET_FILE:Qt5::PrintSupport>
  )
  install(FILES ${QT_DLLS} DESTINATION bin)
elseif(APPLE)
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



install( TARGETS CbctRecon DESTINATION bin )
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/fdk_opencl.cl DESTINATION bin)
# ${CMAKE_CURRENT_BINARY_DIR}/../bin

##-----------------------------------------------------------------------------
##  Additional install files
##-----------------------------------------------------------------------------
if(WIN32)
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

  if(OPENMP_FOUND)
    set(CMAKE_INSTALL_OPENMP_LIBRARIES ON)
  endif(OPENMP_FOUND)
  include(InstallRequiredSystemLibraries)

  if(CUDA_FOUND)
    set (DLL_PATTERNS "cufft64*.dll" "cublas64*.dll" "cudart64*.dll")
    foreach (p ${DLL_PATTERNS})
    file (GLOB DLLFILES "${CUDA_TOOLKIT_ROOT_DIR}/bin/${p}")
    if(DLLFILES)
      install (FILES ${DLLFILES} DESTINATION bin)
    endif()
    endforeach ()
  endif()

  if(ITK_FOUND)
    # ITK*.dll might not include RTK_AS_ITK_EXTERNAL libs
    set (DLL_PATTERNS "*.dll")
    foreach (p ${DLL_PATTERNS})
      file (GLOB DLLFILES "${ITK_DIR}/bin/release/${p}")
      if(DLLFILES)
        install (FILES ${DLLFILES} DESTINATION bin)
      endif()
    endforeach ()
  endif()

  if(RTK_FOUND)
    install (FILES ${RTK_DIR}/bin/release/RTK.dll DESTINATION bin)
    install (FILES ${RTK_DIR}/bin/release/rtkcuda.dll DESTINATION bin)
    install (FILES ${RTK_DIR}/bin/release/ITKCudaCommon.dll DESTINATION bin)
    install (FILES ${RTK_DIR}/bin/release/lpsolve55.dll DESTINATION bin) #is 55 a version number?
  endif()

  if(Plastimatch_FOUND)
    set (DLL_PATTERNS "plm*.dll" "")
    foreach (p ${DLL_PATTERNS})
      file (GLOB DLLFILES "${Plastimatch_BUILD_DIR}/release/${p}")
      if(DLLFILES)
        install (FILES ${DLLFILES} DESTINATION bin)
      endif()
    endforeach ()
  endif()

  if(USE_CLFFT)
    install (FILES ${CLFFT_ROOT_DIR}/bin/clFFT.dll DESTINATION bin)
  endif()

  #Add FFT dlls to Install
  if(FFTW_FOUND)
    if(EXISTS "${FFTW_DIR}/libfftw3-3.dll")
      install (FILES "${FFTW_DIR}/libfftw3-3.dll" DESTINATION bin)
    endif()
    if(EXISTS "${FFTW_DIR}/libfftw3f-3.dll")
      install (FILES "${FFTW_DIR}/libfftw3f-3.dll" DESTINATION bin)
    endif()
  endif()
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
if(CbctRecon_PACKAGE_WIX)
  set (CPACK_GENERATOR "WIX")
elseif(CbctRecon_PACKAGE_NSIS)
  set (CPACK_GENERATOR "NSIS")
else ()
  set (CPACK_GENERATOR "ZIP")
endif()
set (CPACK_SOURCE_GENERATOR "TBZ2")

# For WiX install, and possibly for debian -dev targets,
# we should properly set up the install components.
# Until then, let's just kill it
set (CPACK_MONOLITHIC_INSTALL 1)

# General stuff
set (CPACK_INSTALL_PREFIX "")
set (CPACK_PACKAGE_CONTACT "agravgaard@protonmail.com")
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
set (CPACK_NSIS_HELP_LINK "https://gitlab.com/agravgaard/cbctrecon")
set (CPACK_NSIS_URL_INFO_ABOUT "https://gitlab.com/agravgaard/cbctrecon")
set (CPACK_NSIS_CONTACT "andreasga22@gmail.com")
set (CbctRecon_NSIS_VERSION_STRING "${CbctRecon_VERSION_MAJOR}.${CbctRecon_VERSION_MINOR}.${CbctRecon_VERSION_PATCH}")

# WIX stuff
set (CbctRecon_WINDOWS_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extra/windows-install")
set (CPACK_WIX_HELP_LINK "https://gitlab.com/agravgaard/cbctrecon")
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
if(CPACK_GENERATOR STREQUAL "WIX")
  set (CPACK_PACKAGE_NAME "CbctRecon")
  if(CPACK_PACKAGE_VERSION_TWEAK AND CPACK_PACKAGE_VERSION_TWEAK GREATER "0")
    set (CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.${CPACK_PACKAGE_VERSION_TWEAK}")
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
      set (CPACK_SYSTEM_NAME "win64")
    else ()
      set (CPACK_SYSTEM_NAME "win32")
    endif()
    set (CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME}")
  endif()
  set_property (INSTALL "bin/CbctRecon.exe"
    PROPERTY CPACK_START_MENU_SHORTCUTS "CbctRecon Application")
endif()

# ZIP stuff
if(CPACK_GENERATOR STREQUAL "ZIP")
  set (CPACK_SET_DESTDIR ON)
endif()

include (CPack)

##-----------------------------------------------------------------------------
## Option to generate .clang_complete for hip Vim users using clang_complete
##-----------------------------------------------------------------------------
if (CBCTRECON_CLANG_COMPLETE)
  get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
  set (file_clang_complete "${CMAKE_CURRENT_SOURCE_DIR}/.clang_complete")
  foreach (arg ${inc_dirs})
    set (inc_args "${inc_args} -I${arg}\n")
  endforeach ()
  file (WRITE "${file_clang_complete}" "${inc_args}")
  message (STATUS "Generated ${file_clang_complete}")
endif ()

message(STATUS "Generate vim YCM config file for CbctRecon at ${CMAKE_CURRENT_SOURCE_DIR} ...")
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/.ycm_extra_conf.py.in
    ${CMAKE_CURRENT_SOURCE_DIR}/.ycm_extra_conf.py
)
