include(${CbctRecon_SOURCE_DIR}/cmake/DownloadProject.cmake)

macro(external_dependency NAME URL COMMIT)
    if(${NAME} STREQUAL "Plastimatch")
      set(PATCH_CMD PATCH_COMMAND git apply ${CMAKE_SOURCE_DIR}/External/patches/plm.patch)
    else()
      set(PATCH_CMD PATCH_COMMAND "")
    endif()

    if(NOT TARGET ${NAME})
        message(STATUS "external dependency ${NAME} from ${URL} at ${COMMIT}")
        download_project(
            PROJ "${NAME}"
            GIT_REPOSITORY "${URL}"
            GIT_TAG "${COMMIT}"
            GIT_SHALLOW ON
            ${PATCH_CMD}
        )


    add_subdirectory(${${NAME}_SOURCE_DIR} ${${NAME}_BINARY_DIR})
    else()
        message(STATUS "external dependency ${NAME} already satisfied")
    endif()
endmacro()
