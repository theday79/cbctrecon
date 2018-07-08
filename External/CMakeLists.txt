include(${CbctRecon_SOURCE_DIR}/cmake/externals.cmake)

set(git_protocol https://git)

macro(external_proj NAME)
  if(USE_SYSTEM_${NAME})
    unset(${NAME}_DIR CACHE)
    find_package(${NAME} REQUIRED)
  else()
    include(External_${NAME}.cmake)
    external_dependency(${NAME} ${${NAME}_GIT_REPOSITORY} ${${NAME}_GIT_TAG})
  endif()
endmacro()

# REPOSITORIES
set(DCMTK_GIT_REPOSITORY "${git_protocol}hub.com/DCMTK/DCMTK.git")
set(DCMTK_GIT_TAG master)

set(ITK_GIT_REPOSITORY "${git_protocol}hub.com/InsightSoftwareConsortium/ITK.git")
set(ITK_GIT_TAG "d92873e33e8a54e933e445b92151191f02feab42") #"v4.13.0"

set(plm Plastimatch)
set(${plm}_GIT_REPOSITORY "${git_protocol}lab.com/plastimatch/plastimatch.git")
set(${plm}_GIT_TAG master) #"b9b6c7f59b98d52a93d58288af9d5bc25d00f24b")


# DCMTK
external_proj(DCMTK)
#configure_file(
#  ${CbctRecon_SOURCE_DIR}/External/patches/oflimits.h
#  ${DCMTK_SOURCE_DIR}/ofstd/include/dcmtk/ofstd/
#  COPYONLY
#)

# ITK/RTK
external_proj(ITK)
configure_file(
  ${CbctRecon_SOURCE_DIR}/External/patches/ExternalRTK.remote.cmake
  ${ITK_SOURCE_DIR}/Modules/Remote/
  COPYONLY
)

# Plastimatch
external_proj(${plm})

set(${plm}_DIR ${${plm}_BUILD_DIR})
set(${plm}_SRC ${${plm}_SOURCE_DIR})

#ADD_CUSTOM_COMMAND(
#  TARGET ${proj}
#  PRE_BUILD
#  COMMAND ${CMAKE_COMMAND} -E remove
#    ${proj}/cmake/FindOpenMP.cmake
#  )


