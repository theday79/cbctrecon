include(CTest)
include(ITKExternalData)
include(ITKModuleTest)
include(ITKDownloadSetup)

function(add_cbctrecon_test)
  set(options "")
  set(oneValueArgs
    TNAME
    )
  set(multiValueArgs
    SRC_FILES
    DATA_ARGS
    )
  cmake_parse_arguments(T_ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  add_executable(${T_ARGS_TNAME} ${T_ARGS_SRC_FILES})
  target_link_libraries(${T_ARGS_TNAME}
    PRIVATE CbctReconLib
    )

  target_include_directories(${T_ARGS_TNAME}
    PUBLIC ${CBCTRECON_INCLUDE_DIRS}
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
    )

  find_program(MEMORYCHECK_COMMAND valgrind)

  if(MEMORYCHECK_COMMAND)
    message(STATUS "Running CTest ${T_ARGS_TNAME} with valgrind")
    set(VG_COMM "${MEMORYCHECK_COMMAND} --leak-check=full") 
  endif()
  
  itk_add_test(
    NAME ${T_ARGS_TNAME}
    COMMAND ${VG_COMM} $<TARGET_FILE:${T_ARGS_TNAME}> ${T_ARGS_DATA_ARGS}
    )

endfunction()

enable_testing()
