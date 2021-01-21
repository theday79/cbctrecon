# -----------------------------------------------------------------------------
# HandleITK.cmake Check ITK version and optional components Include use file (for
# registering IO factories)
# -----------------------------------------------------------------------------
# See COPYRIGHT.TXT and LICENSE.TXT for copyright and license information
# -----------------------------------------------------------------------------

# GCS 2016-10-20 The GDCM USE file spews out copious extraneous warnings on Debian. If we
# are not using GDCM, this might be avoided by setting ITK_USE_SYSTEM_GDCM to FALSE.
# However, this has the undesirable effect of removing gdcm includes from include list and
# gdcm libraries from library list. if (DCMTK_FOUND) set (ITK_USE_SYSTEM_GDCM FALSE) endif
# ()

# GCS 2017-12-14 On older ITK version, the use file sets variables such as DCMTK_FOUND,
# DCMTK_DIR.  Needs more investigation.
include(${ITK_USE_FILE})

if(NOT ITK_VERSION)
  set(ITK_VERSION "${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH}")
endif()
if(${ITK_VERSION} VERSION_LESS "3.16.0")
  message(FATAL_ERROR "Fatal Error. ITK must be version 3.16.0 or higher")
endif()
if(${ITK_VERSION_MAJOR} VERSION_EQUAL "3")
  if(NOT ITK_USE_REVIEW)
    message(FATAL_ERROR "Fatal Error. ITK must be compiled with ITK_USE_REVIEW set to ON")
  endif()
  set(ITK_LIBRARIES ${ITK_LIBRARIES} ITKIOReview)
elseif(${ITK_VERSION_MAJOR} VERSION_EQUAL "4")
  if(${ITK_VERSION} VERSION_LESS "4.1")
    message(FATAL_ERROR "Fatal Error. ITK 4 must be 4.1 or greater")
  endif()
else()
  message(AUTHOR_WARNING "Warning. ITK version greater than 4.X is not tested")
endif()
message(STATUS "ITK_VERSION = ${ITK_VERSION} found")

# Find ITK DLL directory.  This is used on Windows for both regression testing and
# packaging.
if(NOT ITK_FOUND)
  set(ITK_BASE "${PLM_BINARY_DIR}/ITK-build")
elseif(${ITK_VERSION} VERSION_LESS "4.1")
  set(ITK_BASE "${ITK_LIBRARY_DIRS}")
else()
  # At some point in time (presumably around ITK 4.1), ITK stopped creating the variable
  # ITK_LIBRARY_DIRS.  Therefore, we infer from the configuration filename. Remove
  # filename
  string(
    REGEX
    REPLACE "/[^/]*$"
            ""
            ITK_LIBRARY_DIRS_41
            "${ITK_CONFIG_TARGETS_FILE}")
  # If configuring against installation directory, walk up to base directory
  string(
    REGEX
    REPLACE "/lib/cmake/ITK-.*]*$"
            ""
            ITK_LIBRARY_DIRS_41
            "${ITK_LIBRARY_DIRS_41}")
  set(ITK_BASE "${ITK_LIBRARY_DIRS_41}")
endif()

message(STATUS "ITK_BASE = ${ITK_BASE}")
if(NOT WIN32)
  set(ITK_DLL_DIR "")
elseif(IS_DIRECTORY "${ITK_BASE}/bin/Release")
  set(ITK_DLL_DIR "${ITK_BASE}/bin/Release")
elseif(IS_DIRECTORY "${ITK_BASE}/Release")
  set(ITK_DLL_DIR "${ITK_BASE}/Release")
elseif(IS_DIRECTORY "${ITK_BASE}/bin")
  set(ITK_DLL_DIR "${ITK_BASE}/bin")
else()
  set(ITK_DLL_DIR "")
endif()
