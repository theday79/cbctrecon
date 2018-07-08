IF (NOT EXISTS "@CMAKE_CURRENT_BINARY_DIR@/install_manifest.txt")
    MESSAGE(FATAL_ERROR "Cannot find install manifest: \"@CMAKE_CURRENT_BINARY_DIR@/install_manifest.txt\"")
ENDIF(NOT EXISTS "@CMAKE_CURRENT_BINARY_DIR@/install_manifest.txt")

file(READ "@CMAKE_CURRENT_BINARY_DIR@/install_manifest.txt" files)
string(REGEX REPLACE "\n" ";" files "${files}")
FOREACH (file ${files})
    MESSAGE(STATUS "Uninstalling \"$ENV{DESTDIR}${file}\"")
    IF (EXISTS "$ENV{DESTDIR}${file}")
        EXECUTE_PROCESS(
            COMMAND @CMAKE_COMMAND@ -E remove "$ENV{DESTDIR}${file}"
            OUTPUT_VARIABLE rm_out
            RESULT_VARIABLE rm_retval
        )
        IF(NOT ${rm_retval} EQUAL 0)
            MESSAGE(FATAL_ERROR "Problem when removing \"$ENV{DESTDIR}${file}\"")
        ENDIF (NOT ${rm_retval} EQUAL 0)
    ELSE (EXISTS "$ENV{DESTDIR}${file}")
        MESSAGE(STATUS "File \"$ENV{DESTDIR}${file}\" does not exist.")
    ENDIF (EXISTS "$ENV{DESTDIR}${file}")
ENDFOREACH(file)
