# Heavily inspired by https://github.com/codeplaysoftware/syclacademy

cmake_minimum_required(VERSION 3.14)

# Options

option(SYCL_USE_COMPUTECPP "Configure to compile with ComputeCpp" OFF)
# option(SYCL_USE_DPCPP "Configure to compile with DPC++" OFF)
option(SYCL_USE_HIPSYCL "Configure to compile with hipSYCL" OFF)
option(SYCL_USE_TRISYCL "Configure to compile with TriSYCL" OFF)

option(USE_SYSTEM_COMPUTECPP_SDK "OFF = Download and configure ComputeCpp SDK"
       ON)

# Find SYCL
include(${CMAKE_CURRENT_SOURCE_DIR}/CMake/SYCLUtils.cmake)

set(CMAKE_FOLDER SYCL)

# Make library

set(syclext_sources syclext_utils.hpp syclext_functions.hpp)

add_library(sycl_ext ${syclext_sources})

add_sycl_to_target(
  TARGET
  sycl_ext
  SOURCES
  ${syclext_sources})
