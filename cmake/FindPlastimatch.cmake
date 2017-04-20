# - Find an Plastimatch installation or build tree.

# When Plastimatch is found, the PlastimatchConfig.cmake file is sourced to setup the
# location and configuration of Plastimatch.  Please read this file, or
# PlastimatchConfig.cmake.in from the Plastimatch source tree for the full list of
# definitions.  Of particular interest is Plastimatch_USE_FILE, a CMake source file
# that can be included to set the include directories, library directories,
# and preprocessor macros.  In addition to the variables read from
# PlastimatchConfig.cmake, this find module also defines
#  Plastimatch_DIR  - The directory containing PlastimatchConfig.cmake.
#             This is either the root of the build tree,
#             or the lib/InsightToolkit directory.
#             This is the only cache entry.
#
#  Plastimatch_FOUND - Whether Plastimatch was found.  If this is true,
#              Plastimatch_DIR is okay.
#

SET(Plastimatch_DIR_STRING "directory containing PlastimatchConfig.cmake.  This is either the root of the build tree, or PREFIX/lib for an installation.")

# Search only if the location is not already known.
IF(NOT Plastimatch_DIR)
  # Get the system search path as a list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" Plastimatch_DIR_SEARCH1 "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" Plastimatch_DIR_SEARCH1 "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/;" ";" Plastimatch_DIR_SEARCH2 ${Plastimatch_DIR_SEARCH1})

  # Construct a set of paths relative to the system search path.
  SET(Plastimatch_DIR_SEARCH "")
  FOREACH(dir ${Plastimatch_DIR_SEARCH2})
    SET(Plastimatch_DIR_SEARCH ${Plastimatch_DIR_SEARCH} "${dir}/../lib")
  ENDFOREACH(dir)

  #
  # Look for an installation or build tree.
  #
  FIND_PATH(Plastimatch_DIR PlastimatchConfig.cmake
    # Look for an environment variable Plastimatch_DIR.
    $ENV{Plastimatch_DIR}

    # Look in places relative to the system executable search path.
    ${Plastimatch_DIR_SEARCH}

    # Look in standard UNIX install locations.
    /usr/local/lib
    /usr/lib

    # Read from the CMakeSetup registry entries.  It is likely that
    # Plastimatch will have been recently built.
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild1]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild2]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild3]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild4]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild5]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild6]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild7]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild8]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild9]
    [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild10]

    # Help the user find it if we cannot.
    DOC "The ${Plastimatch_DIR_STRING}"
  )
ENDIF(NOT Plastimatch_DIR)

# If Plastimatch was found, load the configuration file to get the rest of the
# settings.
IF(Plastimatch_DIR)
  SET(Plastimatch_FOUND 1)
  INCLUDE(${Plastimatch_DIR}/PlastimatchConfig.cmake)
ELSE(Plastimatch_DIR)
  SET(Plastimatch_FOUND 0)
  IF(Plastimatch_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Please set Plastimatch_DIR to the ${Plastimatch_DIR_STRING}")
  ENDIF(Plastimatch_FIND_REQUIRED)
ENDIF(Plastimatch_DIR)
