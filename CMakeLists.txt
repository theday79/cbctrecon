CMAKE_MINIMUM_REQUIRED (VERSION 3.8) # Yes, less won't do.

set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(POLICY CMP0020)
  cmake_policy(SET CMP0020 OLD)
endif()

SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

OPTION(USE_CUDA "Toggle CUDA on or off" OFF) # OFF by default for multiple reasons
OPTION(USE_OPENCL "Toggle OpenCL on or off" ON) # ON by default for multiple reasons

########################START SUPERBUILD STUFF ################################
include(ExternalProject)
include(ExternalProjectDependency)
include(ExternalProjectGenerateProjectDescription)

project(CbctRecon)

FIND_PACKAGE(Git REQUIRED) # We superbuild now.

FIND_PACKAGE(OpenCL REQUIRED)
FIND_PACKAGE(Qt5 COMPONENTS Core Gui OpenGL PrintSupport Widgets REQUIRED)

#-----------------------------------------------------------------------------
# Main application
#-----------------------------------------------------------------------------
if(NOT DEFINED CbctRecon_MAIN_PROJECT)
  set(CbctRecon_MAIN_PROJECT CbctReconApp CACHE INTERNAL "Main project name")
endif()
mark_as_superbuild(CbctRecon_MAIN_PROJECT:STRING)
if(NOT DEFINED ${CbctRecon_MAIN_PROJECT}_APPLICATION_NAME)
  set(${CbctRecon_MAIN_PROJECT}_APPLICATION_NAME CbctRecon CACHE INTERNAL "Main application name")
else()
  if(NOT DEFINED CbctReconApp_APPLICATION_NAME)
    set(CbctReconApp_APPLICATION_NAME CbctRecon)
  endif()
endif()

set(CbctRecon_MAIN_PROJECT_APPLICATION_NAME ${${CbctRecon_MAIN_PROJECT}_APPLICATION_NAME})

if(WIN32)
  option(CbctRecon_BUILD_WIN32_CONSOLE "Build ${PROJECT_NAME} executable as a console app on windows (allows debug output)" ON)
else()
  set(CbctRecon_BUILD_WIN32_CONSOLE OFF)
endif()
mark_as_superbuild(CbctRecon_BUILD_WIN32_CONSOLE:BOOL)
#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(CbctRecon_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif()

foreach(dep QT ${DEPENDENCIES})
  if(CbctRecon_USE_SYSTEM_${dep})
    message(STATUS "Using system ${dep}")
  endif()
endforeach()

set(EXTERNAL_PROJECT_OPTIONAL_ARGS)
if(WIN32)
  list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS -DCbctRecon_SKIP_ROOT_DIR_MAX_LENGTH_CHECK:BOOL=ON)
endif()

#-----------------------------------------------------------------------------
# Superbuild Option - Enabled by default
#-----------------------------------------------------------------------------
option(CbctRecon_SUPERBUILD "Build ${PROJECT_NAME} and the projects it depends on." ON)
mark_as_advanced(CbctRecon_SUPERBUILD)
set(CbctRecon_BINARY_INNER_SUBDIR CbctRecon-build)
#-----------------------------------------------------------------------------
# Sanity checks
#-----------------------------------------------------------------------------
if(WIN32)
  set(${PROJECT_NAME}_ROOT_DIR_MAX_LENGTH 40)
  if(NOT ${PROJECT_NAME}_SUPERBUILD)
    string(LENGTH ${CbctRecon_BINARY_INNER_SUBDIR} _inner_subdir_length)
    math(EXPR ${PROJECT_NAME}_ROOT_DIR_MAX_LENGTH "${${PROJECT_NAME}_ROOT_DIR_MAX_LENGTH} + ${_inner_subdir_length}")
  endif()
  include(PreventDirWithTooManyChars)
endif()
include(PreventInSourceBuilds)
include(PreventInBuildInstalls)
include(PreventDirWithSpaces)
mark_as_superbuild(${PROJECT_NAME}_SKIP_DIR_WITH_SPACES_CHECK:BOOL)
#-----------------------------------------------------------------------------
# Set a default build type if none was specified
#-----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  mark_as_advanced(CMAKE_BUILD_TYPE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()
if(NOT CMAKE_CONFIGURATION_TYPES)
  mark_as_superbuild(VARS CMAKE_BUILD_TYPE ALL_PROJECTS)
endif()

#-----------------------------------------------------------------------------
# Set ITK_REQUIRED_{C,CXX}_FLAGS variables
#-----------------------------------------------------------------------------
include(ITKPlatformSpecificChecks)

#-----------------------------------------------------------------------------
# Set CMAKE_{C,CXX}_FLAGS variables
#-----------------------------------------------------------------------------
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_INIT} ${CbctRecon_REQUIRED_C_FLAGS} ${ITK_REQUIRED_C_FLAGS} ${COVERAGE_C_FLAGS} ${ADDITIONAL_C_FLAGS}" CACHE STRING "CMake C Flags" FORCE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_INIT} ${CbctRecon_REQUIRED_CXX_FLAGS} ${ITK_REQUIRED_CXX_FLAGS} ${COVERAGE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS}" CACHE STRING "CMake CXX Flags" FORCE)

##########################END SUPERBUILD STUFF ################################
