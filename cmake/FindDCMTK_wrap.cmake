##-----------------------------------------------------------------------------
##  See COPYRIGHT.TXT and LICENSE.TXT for copyright and license information
##-----------------------------------------------------------------------------

find_package (DCMTK NO_MODULE REQUIRED)
if (NOT DCMTK_FOUND)
  message (STATUS "Searching for DCMTK using legacy method")
  find_package (DCMTK_legacy REQUIRED)
endif ()

# The DCMTK 3.6.2 DCMTKConfig.cmake seems to be broken on windows
string (REPLACE "DCMTK_INCLUDE_DIRS-NOTFOUND;" "" DCMTK_INCLUDE_DIRS
  "${DCMTK_INCLUDE_DIRS}")
