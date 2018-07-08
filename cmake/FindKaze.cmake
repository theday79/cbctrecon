# - Find libkaze
#  KAZE_INCLUDE_DIR - where to find kaze.h, etc.
#  KAZE_LIBRARIES   - List of libraries when using kaze.
#  KAZE_FOUND       - True if kaze found.

set (KAZE_DIR "" CACHE PATH "Root of libkaze install tree (optional).")

if (KAZE_INCLUDE_DIR)
  # Already in cache, be silent
  set (Kaze_FIND_QUIETLY TRUE)
endif ()


find_path (KAZE_INCLUDE_DIR kimage.h
  ${KAZE_DIR}/include)

set (KAZE_NAMES kaze)
find_library (KAZE_LIBRARY NAMES ${KAZE_NAMES}
  PATHS
  ${KAZE_DIR}/lib)

# handle the QUIETLY and REQUIRED arguments and set KAZE_FOUND to TRUE if 
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Kaze DEFAULT_MSG 
  KAZE_LIBRARY 
  KAZE_INCLUDE_DIR)

if (KAZE_FOUND)
  set (KAZE_LIBRARIES ${KAZE_LIBRARY})
else ()
  set (KAZE_LIBRARIES)
endif ()

mark_as_advanced (KAZE_LIBRARY KAZE_INCLUDE_DIR)
