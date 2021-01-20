# JAS 08.19.2010 I have tested this working under Linux using gcc-4.4 The OSX/Darwin case
# is untested, but should work in theory.

# LINUX: We check for SSE extensions using a RegEx on /proc/cpuinfo
if(CMAKE_SYSTEM_NAME MATCHES "Linux")

  # Store /proc/cpuinfo output into CPUINFO
  exec_program(
    cat ARGS
    "/proc/cpuinfo"
    OUTPUT_VARIABLE CPUINFO)

  # Check for SSE2
  string(
    REGEX
    REPLACE "^.*(sse2).*$"
            "\\1"
            SSE_THERE
            ${CPUINFO})
  string(
    COMPARE EQUAL
            "sse2"
            "${SSE_THERE}"
            SSE2_TRUE)
  if(SSE2_TRUE)
    set(SSE2_FOUND
        true
        CACHE BOOL "SSE2 Available?")
    set(SSE2_FLAGS "-msse2 -mfpmath=sse")
  else(SSE2_TRUE)
    set(SSE2_FOUND
        false
        CACHE BOOL "SSE2 Available?")
  endif(SSE2_TRUE)

  # OSX/DARWIN: We check for SSE extensions using a RegEx on sysctl output
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  exec_program("/usr/sbin/sysctl -n machdep.cpu.features" OUTPUT_VARIABLE CPUINFO)

  # Check for SSE2
  string(
    REGEX
    REPLACE "^.*(SSE2).*$"
            "\\1"
            SSE_THERE
            ${CPUINFO})
  string(
    COMPARE EQUAL
            "SSE2"
            "${SSE_THERE}"
            SSE2_TRUE)
  if(SSE2_TRUE)
    set(SSE2_FOUND
        true
        CACHE BOOL "SSE2 Available?")
    set(SSE2_FLAGS "-msse2 -mfpmath=sse") # Darwin uses gcc, right?
  else(SSE2_TRUE)
    set(SSE2_FOUND
        false
        CACHE BOOL "SSE2 Available?")
  endif(SSE2_TRUE)

  # WINDOWS: Currently, no SSE detection method.  Disabled.
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
  set(SSE2_FOUND
      false
      CACHE BOOL "SSE2 Available?")

  # Something else... BSD or BeOS, perhaps?  Disabled.
else(CMAKE_SYSTEM_NAME MATCHES "Linux")
  set(SSE2_FOUND
      false
      CACHE BOOL "SSE2 Available?")

endif(CMAKE_SYSTEM_NAME MATCHES "Linux")

# Report SSE2 compiler flags or failure.
if(SSE2_FOUND)
  message(STATUS "SSE2_FLAGS \"${SSE2_FLAGS}\"")
else(SSE2_FOUND)
  message(STATUS "CPU does not support SSE2.")
endif(SSE2_FOUND)

# Put this in the advanced toggles
mark_as_advanced(SSE2_FOUND)
