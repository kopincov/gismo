 ### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

PROJECT(QtGismo)

SET( ${PROJECT_NAME}_SOURCES  main.cpp gsMainWindow.cpp )
SET( ${PROJECT_NAME}_HEADERS gsMainWindow.h)
SET( ${PROJECT_NAME}_FORMS gsMainWindow.ui)
#SET( ${PROJECT_NAME}_RESOURCES images.qrc)


# enable warnings
ADD_DEFINITIONS( -Wall )


# by default only QtCore and QtGui modules are enabled
# other modules must be enabled like this:
#SET( QT_USE_QT3SUPPORT TRUE )   
#SET( QT_USE_QTXML TRUE )

# this command finds Qt4 libraries and sets all required variables
# note that it's Qt4, not QT4 or qt4

FIND_PACKAGE( Qt4 REQUIRED )


# add some useful macros and variables
# (QT_USE_FILE is a variable defined by FIND_PACKAGE( Qt4 ) that contains a path to CMake script)
INCLUDE( ${QT_USE_FILE} )



# this command will generate rules that will run rcc on all files from SAMPLE_RCS
# in result SAMPLE_RC_SRCS variable will contain paths to files produced by rcc
#QT4_ADD_RESOURCES( SAMPLE_RC_SRCS ${SAMPLE_RCS} )


# this will run uic on .ui files:
QT4_WRAP_UI(${PROJECT_NAME}_FORMS_HEADERS ${${PROJECT_NAME}_FORMS})

# and finally this will run moc:
QT4_WRAP_CPP(${PROJECT_NAME}_HEADERS_MOC ${${PROJECT_NAME}_HEADERS})

# we need this to be able to include headers produced by uic in our code
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

#QT4_ADD_RESOURCES(helloworld_RESOURCES_RCC ${helloworld_RESOURCES})

ADD_DEFINITIONS(${QT_DEFINITIONS})


ADD_EXECUTABLE(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES}
    ${${PROJECT_NAME}_HEADERS}	
    ${${PROJECT_NAME}_HEADERS_MOC}
    ${${PROJECT_NAME}_FORMS_HEADERS}
)

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${QT_LIBRARIES})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
