######################################################
##  Check if nullptr is supported by the compiler
##  Returns: HAVE_NULLPTR
######################################################
include (CheckCXXSourceCompiles)
set (NULLPTR_TEST_SOURCE "int main() {(void)nullptr;}")
check_cxx_source_compiles ("${NULLPTR_TEST_SOURCE}" HAVE_NULLPTR)
if (NOT HAVE_NULLPTR)
  unset (HAVE_NULLPTR CACHE)
  push_var (CMAKE_REQUIRED_FLAGS)
  set (CMAKE_REQUIRED_FLAGS "-std=c++11")
  check_cxx_source_compiles ("${NULLPTR_TEST_SOURCE}" HAVE_NULLPTR)
  pop_var (CMAKE_REQUIRED_FLAGS)
  if (HAVE_NULLPTR)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  endif ()
endif ()
