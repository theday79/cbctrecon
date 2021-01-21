# ########################################################################################
# OpenMP
# ########################################################################################
include(CheckFunctionExists)
message(STATUS "Check for compiler OpenMP support...")
set(OPENMP_FLAGS)
set(OPENMP_LIBRARIES)
set(OPENMP_FOUND FALSE)

# Key: CFLAGS##LDFLAGS#LIBRARIES Neither CFLAGS nor LDFLAGS can be empty.  Use NONE
# instead.
set(OPENMP_FLAGS_AND_LIBRARIES
    # gcc
    "-fopenmp##-fopenmp#"
    "-fopenmp##-fopenmp#gomp"
    "-fopenmp##-fopenmp#gomp pthread"
    # icc
    "-openmp##-openmp#"
    "-openmp -parallel##-openmp -parallel#"
    # SGI & PGI
    "-mp##-mp#"
    # Sun
    "-xopenmp##-xopenmp#"
    # Tru64
    "-omp##-omp#"
    # AIX
    "-qsmp=omp##-qsmp=omp#"
    # MSVC
    "/openmp##NONE#")

# Massive hack to workaround CMake limitations
list(LENGTH OPENMP_FLAGS_AND_LIBRARIES NUM_FLAGS)
math(EXPR NUM_FLAGS "${NUM_FLAGS} - 1")
foreach(I RANGE 0 ${NUM_FLAGS})
  if(NOT OPENMP_FOUND)
    list(
      GET
      OPENMP_FLAGS_AND_LIBRARIES
      ${I}
      TMP)
    string(
      REGEX MATCH
            "([^#]*)"
            OPENMP_FLAGS
            ${TMP})
    string(
      REGEX
      REPLACE "[^#]*##"
              ""
              TMP
              ${TMP})
    string(
      REGEX MATCH
            "([^#]*)"
            OPENMP_LDFLAGS
            ${TMP})
    string(
      REGEX
      REPLACE "[^#]*#"
              ""
              OPENMP_LIBRARIES
              ${TMP})
    # MESSAGE(STATUS "OPENMP_FLAGS=${OPENMP_FLAGS}") MESSAGE(STATUS "OPENMP_LDFLAGS =
    # ${OPENMP_LDFLAGS}") MESSAGE(STATUS "OPENMP_LIBRARIES = ${OPENMP_LIBRARIES}")
    # MESSAGE(STATUS "-------")

    if(OPENMP_LDFLAGS MATCHES "NONE")
      set(OPENMP_LDFLAGS "")
    endif(OPENMP_LDFLAGS MATCHES "NONE")
    if(OPENMP_LIBRARIES MATCHES " ")
      string(
        REPLACE " "
                ";"
                OPENMP_LIBRARIES
                ${OPENMP_LIBRARIES})
    endif(OPENMP_LIBRARIES MATCHES " ")

    # I think I need to do a try-compile
    set(CMAKE_REQUIRED_FLAGS ${OPENMP_FLAGS})
    set(CMAKE_REQUIRED_LIBRARIES ${OPENMP_LIBRARIES})
    check_function_exists(omp_get_thread_num OPENMP_FOUND${I})

    if(OPENMP_FOUND${I})
      set(OPENMP_FOUND TRUE)
    endif(OPENMP_FOUND${I})
  endif(NOT OPENMP_FOUND)
endforeach(
  I
  RANGE
  0
  ${NUM_FLAGS})

if(OPENMP_FOUND)
  message(
    STATUS "OpenMP flags \"${OPENMP_FLAGS}\", OpenMP libraries \"${OPENMP_LIBRARIES}\"")
else(OPENMP_FOUND)
  message(STATUS "Given compiler does not support OpenMP.")
endif(OPENMP_FOUND)
