##-----------------------------------------------------------------------------
##  HandleITK.cmake
##    Check ITK version and optional components
##    Include use file (for registering IO factories)
##-----------------------------------------------------------------------------
##  See COPYRIGHT.TXT and LICENSE.TXT for copyright and license information
##-----------------------------------------------------------------------------

# GCS 2016-10-20
# The GDCM USE file spews out copious extraneous warnings on Debian.
# If we are not using GDCM, this might be avoided by setting
# ITK_USE_SYSTEM_GDCM to FALSE.  However, this has the undesirable
# effect of removing gdcm includes from include list and gdcm libraries
# from library list.
#if (DCMTK_FOUND)
# set (ITK_USE_SYSTEM_GDCM FALSE)
#endif ()

# GCS 2017-12-14 On older ITK version, the use file sets variables such
# as DCMTK_FOUND, DCMTK_DIR.  Needs more investigation.
include (${ITK_USE_FILE})

if (NOT ITK_VERSION)
  set (ITK_VERSION
    "${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH}")
endif ()
if (${ITK_VERSION} VERSION_LESS "3.16.0")
  message (FATAL_ERROR 
    "Fatal Error. ITK must be version 3.16.0 or higher")
endif ()
if (${ITK_VERSION_MAJOR} VERSION_EQUAL "3")
  if (NOT ITK_USE_REVIEW)
    message (FATAL_ERROR 
      "Fatal Error. ITK must be compiled with ITK_USE_REVIEW set to ON")
  endif ()
  set (ITK_LIBRARIES ${ITK_LIBRARIES} ITKIOReview)
elseif (${ITK_VERSION_MAJOR} VERSION_EQUAL "4")
  if (${ITK_VERSION} VERSION_LESS "4.1")
    message (FATAL_ERROR 
      "Fatal Error. ITK 4 must be 4.1 or greater")
  endif ()
else ()
  message (STATUS 
    "What if ITK not version 3.X or 4.X?")
endif ()
message (STATUS "ITK_VERSION = ${ITK_VERSION} found")
