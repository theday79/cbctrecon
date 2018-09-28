include(CTest)
include(${ITK_SOURCE_DIR}/CMake/ITKExternalData.cmake)
include(${ITK_SOURCE_DIR}/CMake/ITKModuleTest.cmake)
include(${ITK_SOURCE_DIR}/CMake/ITKDownloadSetup.cmake)

function(add_cbctrecon_test)
  cmake_parse_arguments(
    ARGS
    ""
    "TARGET"
    "SRC_FILES;DATA_ARGS"
    ${ARGN}
    )

  add_executable(${ARGS_TARGET} ${ARGS_SRC_FILES})
  target_link_libraries(${ARGS_TARGET}
    PRIVATE CbctReconLib
    )

  target_include_directories(${ARGS_TARGET}
    PUBLIC ${CBCTRECON_INCLUDE_DIRS}
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
    )

  find_program(MEMORYCHECK_COMMAND valgrind)

  if(MEMORYCHECK_COMMAND)
    message(STATUS "Running CTest ${ARGS_TARGET} with valgrind")
    set(VG_COMM ${MEMORYCHECK_COMMAND} "--leak-check=full") 
  endif()
  
  itk_add_test(
    NAME ${ARGS_TARGET}
    COMMAND ${VG_COMM} $<TARGET_FILE:${ARGS_TARGET}> ${ARGS_DATA_ARGS}
    )

endfunction()

enable_testing()
