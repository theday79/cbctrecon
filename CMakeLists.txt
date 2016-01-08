CMAKE_MINIMUM_REQUIRED (VERSION 2.8)


# This project is designed to be built outside the RTK source tree.


#test
PROJECT(CbctRecon)


# Find the RTK libraries and includes


SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

FIND_PACKAGE(RTK REQUIRED)
INCLUDE(${RTK_USE_FILE})

FIND_PACKAGE(CUDA)


get_filename_component(MY_PARENT_DIR ${CUDA_LIBRARIES} PATH)
set (CUDA_FFT_LIB ${MY_PARENT_DIR}/cufft.lib)

FIND_PACKAGE(Qt4 REQUIRED)

find_package (Plastimatch)
message (STATUS "Plastimatch_FOUND = ${Plastimatch_FOUND}")
message (STATUS "Plastimatch_DIR = ${Plastimatch_DIR}")
message (STATUS "Plastimatch_SOURCE = ${Plastimatch_SRC}")
message (STATUS "PLM_SRC_ROOT = ${PLM_SRC_ROOT}")

message (STATUS "CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}")
message (STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}")
message (STATUS "MY_PARENT_DIR = ${MY_PARENT_DIR}")

message (STATUS "PLASTIMATCH_LIBRARIES = ${PLASTIMATCH_LIBRARIES}")
message (STATUS "PLASTIMATCH_LIBRARY = ${PLASTIMATCH_LIBRARY}")
message (STATUS "PLASTIMATCH_INCLUDE_DIRS = ${PLASTIMATCH_INCLUDE_DIRS}")
message (STATUS "PLASTIMATCH_INCLUDE_DIRS = ${PLASTIMATCH_INCLUDE_DIRS}")
message (STATUS "PLASTIMATCH_SOURCE_DIR = ${PLASTIMATCH_DIR}")

#MESSAGE(STATUS "PLASTIMATCH_LINK_DIRECTORIES = ${Plastimatch_LINK_DIRECTORIES}")
#MESSAGE(STATUS "PLASTIMATCH_LINK_DIRECTORIES = ${Plastimatch_LINK_DIRS}")
#MESSAGE(STATUS "PLASTIMATCH_SOURCE_DIRECTORY = ${PLASTIMATCH_LINK_DIR}")

if (QT4_FOUND)
  message (STATUS "Looking for Qt4 - found")
  include (${QT_USE_FILE})
else ()
  message (STATUS "Looking for Qt4 - not found")
endif ()

#Find ITk folder
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
   INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR 
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

#Find DCMTK

find_package(DCMTK REQUIRED)
#LINK_DIRECTORIES(${DCMTK_LIBRARY_DIRS}) #not defined
#LINK_DIRECTORIES("H:/lib/dcmtk3.6.1/BUILD32/lib/Debug") #works
LINK_DIRECTORIES("${DCMTK_DIR}/lib")

#works But A lot of lInk error (detection...

MESSAGE(STATUS "DCMTK_FOUND = ${DCMTK_FOUND}")
MESSAGE(STATUS "DCMTK_DIR = ${DCMTK_DIR}")
MESSAGE(STATUS "DCMTK_ROOT_DIR = ${DCMTK_ROOT_DIR}")
MESSAGE(STATUS "DCMTK_INCLUDE_DIRS = ${DCMTK_INCLUDE_DIRS}")
MESSAGE(STATUS "DCMTK_LINK_DIRECTORIES = ${DCMTK_LINK_DIRECTORIES}")
MESSAGE(STATUS "DCMTK_LIBRARIES = ${DCMTK_LIBRARIES}")

message (STATUS "DCMTK_SOURCE_DIRECTORY = ${DCMTK_SOURCE_DIRECTORY}")


SET(CbctRecon_SOURCES
	main.cpp
	cbctrecon.cpp
	YK16GrayImage.cxx
	qyklabel.cpp
	qcustomplot.cpp
	DlgRegistration.cxx
	DlgExternalCommand.cxx
	)

SET(CbctRecon_HEADERS
	cbctrecon.h
	YK16GrayImage.h
	qyklabel.h
	qcustomplot.h
	DlgRegistration.h
	DlgExternalCommand.h
	)

SET(CbctRecon_FORMS
	cbctrecon.ui
	DlgRegistration.ui
	DlgExternalCommand.ui)



INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

QT4_WRAP_CPP(CbctRecon_HEADERS_MOC ${CbctRecon_HEADERS})
QT4_WRAP_UI(CbctRecon_FORMS_HEADERS ${CbctRecon_FORMS})
QT4_ADD_RESOURCES(CbctRecon_RESOURCES_RCC ${CbctRecon_RESOURCES})


# Executable


ADD_EXECUTABLE(CbctRecon ${CbctRecon_SOURCES} 
	${CbctRecon_HEADERS} 
    ${CbctRecon_HEADERS_MOC} 
    ${CbctRecon_FORMS_HEADERS} 
     )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

#SET (PLASTIMATCH_INCLUDE_DIRS
#	${PLASTIMATCH_INCLUDE_DIRS}
#	"H:/lib/plastimatch/SRC/src/plastimatch/segment"
#	"H:/lib/plastimatch/SRC/src/plastimatch/cli"
#	"H:/lib/plastimatch/SRC/src/plastimatch/clp")

#SET (PLASTIMATCH_INCLUDE_DIRS
#	${PLASTIMATCH_INCLUDE_DIRS})

message (STATUS "PLASTIMATCH_SOURCE_DIRECTORY = ${PLASTIMATCH_INCLUDE_DIRS}")
message (STATUS "CMAKE_LIBRARY_PATH = ${CMAKE_LIBRARY_PATH}")
message (STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")
	

#this is for "C++ General-> Add below to Additional Include Directories"
INCLUDE_DIRECTORIES(${PLASTIMATCH_INCLUDE_DIRS})

SET(QT_USE_QTCORE TRUE)
SET(QT_USE_QTGUI TRUE)
SET(QT_USE_QTOPENGL TRUE)

#This is for "Linker->additional library directories" for DCMRT dir (to deal with oflog.lib not found error)
#LINK_DIRECTORIES("H:/lib/dcmtk3.6.1/BUILD32/lib/Debug")
#SET_TARGET_PROPERTIES (CbctRecon PROPERTIES LINK_FLAGS "H:\lib\dcmtk3.6.1\BUILD32\lib\Debug")

message (STATUS "RTK_USE_FILE = ${RTK_USE_FILE}")

message (STATUS "RTK_LIBRARY_DIRS = ${RTK_LIBRARY_DIRS}")
message (STATUS "DCMTK_LIBRARY_DIRS = ${DCMTK_LIBRARY_DIRS}")
message (STATUS "DCMTK_LIBRARIES = ${DCMTK_LIBRARIES}")

#this is for "Linker->Additional dependencies"
TARGET_LINK_LIBRARIES(CbctRecon
netapi32.lib
${ITK_LIBRARIES}
${RTK_LIBRARIES}
${QT_LIBRARIES}
${QT_QTCORE_LIBRARY}
${QT_QTGUI_LIBRARY}
${QT_QTOPENGL_LIBRARY}
${CUDA_LIBRARIES}
${CUDA_FFT_LIB}
${PLASTIMATCH_LIBRARIES}
#${DCMTK_LIBRARIES} a lot of link error.. dectection.. 
)

