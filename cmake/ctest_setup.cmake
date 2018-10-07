include(CTest)
get_filename_component(_CBCTRECONExternalData_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(ExternalData)
set(ExternalData_URL_ALGO_MD5_lower md5)
set(ExternalData_URL_TEMPLATES
  # Data published on Girder
  "https://data.kitware.com/api/v1/file/hashsum/%(algo)/%(hash)/download"
  )

set(ExternalData_LINK_CONTENT MD5)

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
  
  ExternalData_Add_test(CbctData
    NAME ${ARGS_TARGET}
    COMMAND ${VG_COMM} $<TARGET_FILE:${ARGS_TARGET}> ${ARGS_DATA_ARGS}
    )

endfunction()

enable_testing()
