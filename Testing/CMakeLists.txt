set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

function(unpack_ext_data)
  cmake_parse_arguments(
    ARGS
    ""
    "TARGET;TAR_FILE;TAR_OUTPUT_PATH;TAR_SUBOUTPUT"
    ""
    ${ARGN}
    )

  add_custom_command(
    TARGET ${ARGS_TARGET}
    POST_BUILD # Because after build the tars must have been downloaded
    COMMAND ${CMAKE_COMMAND} -E tar xzf ${ARGS_TAR_FILE}
    DEPENDS ${ARGS_TAR_FILE}
    WORKING_DIRECTORY ${ARGS_TAR_OUTPUT_PATH}
    COMMENT "Unpacking ${ARGS_TAR_FILE}"
    VERBATIM
  )
endfunction()


set(test_target
  cbctrecon_test
  )

set(${test_target}_SRC
  cbctrecon_test.cpp
  cbctrecon_test.hpp
  )

add_library(${test_target}
  ${${test_target}_SRC})

target_link_libraries(${test_target}
  PUBLIC CbctReconLib
  )

target_include_directories(${test_target}
  PUBLIC ${CBCTRECON_INCLUDE_DIRS}
  PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
  )


if(CMAKE_COMPILER_IS_GNUCXX AND CBCTRECON_COVERAGE)
  include(CodeCoverage)
  setup_target_for_coverage(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    NAME CbctReconLib_coverage
    SOURCES Library/CbctReconLib/*
    COLLECT_EXISTING
    )
  target_link_libraries(${test_target} PRIVATE gcov)
endif()

# Test: Reading Dicom files

set(test_target
  test_dcm_reader
  )

set(${test_target}_DATA
  DATA{Input/Alderson/CT.tar.gz}
  )

add_cbctrecon_test(
  TARGET ${test_target}
  DATA_ARGS ${${test_target}_DATA}
  )

set(TAR_OUTPUT_PATH ${CMAKE_BINARY_DIR}/ExternalData/Testing/Input/Alderson)
set(TAR_FILE ${TAR_OUTPUT_PATH}/CT.tar.gz)
set(TAR_SUBOUTPUT CT)

unpack_ext_data(
  TARGET ${test_target}
  TAR_FILE ${TAR_FILE}
  TAR_OUTPUT_PATH ${TAR_OUTPUT_PATH}
  TAR_SUBOUTPUT ${TAR_SUBOUTPUT}
  )

# Test: CB projection reader

set(test_target
    test_cb_reader
    )

set(${test_target}_DATA
  DATA{Input/Alderson/CBCT.tar.gz}
  DATA{Baseline/Alderson/CBproj.mha}
  )

add_cbctrecon_test(
  TARGET ${test_target}
  DATA_ARGS ${${test_target}_DATA}
  )

set(TAR_FILE ${TAR_OUTPUT_PATH}/CBCT.tar.gz)
set(TAR_SUBOUTPUT CBCT)

unpack_ext_data(
  TARGET ${test_target}
  TAR_FILE ${TAR_FILE}
  TAR_OUTPUT_PATH ${TAR_OUTPUT_PATH}
  TAR_SUBOUTPUT ${TAR_SUBOUTPUT}
  )

# Test: End 2 End

set(test_target
    test_end_to_end
    )

set(${test_target}_DATA
  DATA{Input/Alderson/CT.tar.gz}
  DATA{Input/Alderson/CBCT.tar.gz}
  )
  
add_cbctrecon_test(
  TARGET ${test_target}
  DATA_ARGS ${${test_target}_DATA}
  )
