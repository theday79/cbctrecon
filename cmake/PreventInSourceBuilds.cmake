# This function will prevent in-source builds. Copied from ITK.
function(AssureOutOfSourceBuilds)
  # make sure the user doesn't play dirty with symlinks
  get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

  # disallow in-source builds
  if("${srcdir}" STREQUAL "${bindir}")
    message(FATAL_ERROR "Error: Plastimatch should not be built in its source directory")
  endif()
endfunction()

AssureOutOfSourceBuilds()
