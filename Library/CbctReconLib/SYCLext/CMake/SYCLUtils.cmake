# Variables

set(SYCL_INSTALL_ROOT CACHE STRING "NOT-FOUND")

set(SYCL_IMPLEMENTATIONS_USED 0)
if(SYCL_USE_COMPUTECPP)
  math(EXPR SYCL_IMPLEMENTATIONS_USED ${SYCL_IMPLEMENTATIONS_USED}+1)
endif()
if(SYCL_USE_DPCPP)
  math(EXPR SYCL_IMPLEMENTATIONS_USED ${SYCL_IMPLEMENTATIONS_USED}+1)
endif()
if(SYCL_USE_HIPSYCL)
  math(EXPR SYCL_IMPLEMENTATIONS_USED ${SYCL_IMPLEMENTATIONS_USED}+1)
endif()
if(SYCL_USE_TRISYCL)
  math(EXPR SYCL_IMPLEMENTATIONS_USED ${SYCL_IMPLEMENTATIONS_USED}+1)
endif()
if(${SYCL_IMPLEMENTATIONS_USED} EQUAL 0)
  message(
    FATAL_ERROR
      "No SYCL implementation specified, please set either "
      "SYCL_USE_COMPUTECPP or SYCL_USE_HIPSYCL or SYCL_USE_TRISYCL to ON.")
endif()
if(${SYCL_IMPLEMENTATIONS_USED} GREATER 1)
  message(
    FATAL_ERROR
      "Multiple SYCL implementations specified, please only "
      "set one of SYCL_USE_COMPUTECPP or  SYCL_USE_HIPSYCL or SYCL_USE_TRISYCL to ON."
  )
endif()

# Common setup

include(${CbctRecon_SOURCE_DIR}/cmake/externals.cmake)
set(git_protocol https://git)

macro(external_proj NAME)
  if(USE_SYSTEM_${NAME})
    find_package(${NAME} QUIET)
  else()
    include(External_${NAME}.cmake)
    external_dependency(${NAME} ${${NAME}_GIT_REPOSITORY} ${${NAME}_GIT_TAG})
  endif()
endmacro()

# REPOSITORIES
set(COMPUTECPP_SDK_GIT_REPOSITORY
    "${git_protocol}hub.com/codeplaysoftware/computecpp-sdk.git")
set(COMPUTECPP_SDK_GIT_TAG master)

external_proj(COMPUTECPP_SDK)

list(APPEND CMAKE_MODULE_PATH "${COMPUTECPP_SDK_SOURCE_DIR}/cmake/Modules/")

# ComputeCpp setup

if(SYCL_USE_COMPUTECPP)
  if(NOT SYCL_INSTALL_ROOT)
    message(
      FATAL_ERROR
        "SYCL implementation root not provided, please specify "
        "the path to the root of the chosen SYCL implementation using "
        "SYCL_INSTALL_ROOT=<path/to/install/root>.")
  endif()

  set(ComputeCpp_DIR ${SYCL_INSTALL_ROOT})
  include(External/computecpp-sdk/cmake/Modules/ComputeCppCompilerChecks.cmake)
  find_package(ComputeCpp REQUIRED)
endif()

# DPC++ setup

if(SYCL_USE_DPCPP)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    check_cxx_compiler_flag("-fsycl" sycl_support)
    if(sycl_support)
      message(STATUS "Using clang as SYCL compiler")
    else()
      message(FATAL_ERROR "This version of clang does not support SYCL")
    endif()
  else()
    message(FATAL_ERROR "CXX compileri ID must be Clang for dpcpp")
  endif()

  enable_language(SYCL)

  function(add_sycl_to_target)
    set(options)
    set(one_value_keywords TARGET)
    set(multi_value_keywords SOURCES)
    cmake_parse_arguments(ADD_SYCL "${options}" "${one_value_keywords}"
                          "${multi_value_keywords}" ${ARGN})

    if(NOT OpenCL_FOUND)
      find_package(OpenCL REQUIRED)
    endif()
    target_link_libraries(${ADD_SYCL_TARGET} PUBLIC OpenCL::OpenCL)
    target_compile_options(${ADD_SYCL_TARGET} PUBLIC "-fsycl")

    set_source_files_properties(${ADD_SYCL_SOURCES} PROPERTIES LANGUAGE SYCL)
  endfunction()

endif()

# hipSYCL setup

if(SYCL_USE_HIPSYCL)
  find_package(hipSYCL CONFIG REQUIRED PATHS ${SYCL_INSTALL_ROOT}/lib/cmake
               /opt/hipSYCL/ROCm/lib/cmake)
endif()

# TriSYCL setup

if(SYCL_USE_TRISYCL)
  find_package(TriSYCL CONFIG REQUIRED PATHS ${SYCL_INSTALL_ROOT}/lib/cmake)
endif()
