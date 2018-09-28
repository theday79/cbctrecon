include(ExternalData)

function(add_cbctrecon_test NAME)
    find_program(MEMORYCHECK_COMMAND valgrind)

    if(MEMORYCHECK_COMMAND)
        message(STATUS "Running CTest ${NAME} with valgrind")
	ExternalData_add_test(${NAME}_Data
            NAME ${NAME}
	    COMMAND ${MEMORYCHECK_COMMAND} --leak-check=full $<TARGET_FILE:${NAME}> ${ARGN})
    else()
	    ExternalData_add_test(CbctReconData
                NAME ${NAME}
		COMMAND $<TARGET_FILE:${NAME}> ${ARGN})
    endif()
endfunction()

enable_testing()
