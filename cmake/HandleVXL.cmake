##-----------------------------------------------------------------------------
##  HandleVXL.cmake
##    Ensure that VXL is working
##-----------------------------------------------------------------------------
##  See COPYRIGHT.TXT and LICENSE.TXT for copyright and license information
##-----------------------------------------------------------------------------
include (CheckCXXSourceCompiles)

set (VXL_TEST_SOURCE 
  "#include <itkImage.h>
   int main (int argc, char* argv[]) {return 0;}")

push_vars ("CMAKE_REQUIRED_DEFINITIONS" "CMAKE_REQUIRED_INCLUDES")
set (CMAKE_REQUIRED_INCLUDES ${ITK_INCLUDE_DIRS})
    set (CMAKE_REQUIRED_QUIET TRUE)
check_cxx_source_compiles ("${VXL_TEST_SOURCE}" VXL_IS_OK)
if (NOT VXL_IS_OK)
  set (CMAKE_REQUIRED_DEFINITIONS "-DVCL_CAN_STATIC_CONST_INIT_FLOAT=0")
  check_cxx_source_compiles ("${VXL_TEST_SOURCE}" VXL_NEEDS_FIXING)
  if (VXL_NEEDS_FIXING)
    add_definitions ("-DVCL_CAN_STATIC_CONST_INIT_FLOAT=0")
  endif ()
endif ()
pop_vars ("CMAKE_REQUIRED_DEFINITIONS" "CMAKE_REQUIRED_INCLUDES")

