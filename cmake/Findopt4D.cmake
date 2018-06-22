# - Find opt4D
#
#  opt4D_INCLUDE_DIR
#  opt4D_LIBRARIES
#  opt4D_FOUND


if (opt4D_INCLUDE_DIR)
  # Already in cache, be silent
  set (opt4D_FIND_QUIETLY TRUE)
endif ()

if (opt4D_DIR)
  find_library (opt4D_LIBRARIES opt4D PATHS ${opt4D_DIR})
endif ()
if (opt4D_LIBRARIES)
  set (opt4D_cache_file "${opt4D_DIR}/CMakeCache.txt")
  if (EXISTS "${opt4D_cache_file}")
    file (READ "${opt4D_cache_file}" opt4D_cache_contents)
    string (REGEX MATCH "opt4D_SOURCE_DIR[^=]*=([^\n]*)\n" opt4D_SOURCE_DIR
      "${opt4D_cache_contents}")
    find_path (opt4D_INCLUDE_DIR BixelVector.hpp PATHS "${CMAKE_MATCH_1}/src"
      NO_DEFAULT_PATH)
  endif ()
endif ()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (opt4D DEFAULT_MSG 
  opt4D_LIBRARIES
  opt4D_INCLUDE_DIR)

mark_as_advanced (opt4D_LIBRARIES opt4D_INCLUDE_DIR)
