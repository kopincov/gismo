######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012-2015 - RICAM-Linz.
######################################################################

cmake_minimum_required(VERSION 2.8.8)

IF( COMMAND CMAKE_POLICY )
  CMAKE_POLICY( SET CMP0003 NEW )
ENDIF( COMMAND CMAKE_POLICY )

# Info:
# http://www.paraview.org/Wiki/ParaView/Users_Guide/Loading_Plugins

project(gsParaview)

IF (ParaView_SOURCE_DIR)
  INCLUDE_DIRECTORIES(
    ${PARAVIEW_INCLUDE_DIRS}
    ${PARAVIEW_GUI_INCLUDE_DIRS}
    ${PARAVIEW_KWSYS_INCLUDE_DIRS}
  )
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF (ParaView_SOURCE_DIR)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/plParaview)

ADD_PARAVIEW_PLUGIN(GiSmoPlugin "1.0"
  SERVER_MANAGER_XML gsParaview.xml
  SERVER_MANAGER_SOURCES vtkGismoReader.cxx
  #GUI_RESOURCE_FILES GismoReaderGUI.xml # deprecated
  )

#If shared library and cpp files #include the headers of another
#library (Say, QtNetwork for example), but your header files don't
#include QtNetwork headers, then QtNetwork is a PRIVATE dependency.
#If your source files and your headers include the headers of another
#library, then it is a PUBLIC dependency.
#If your header files but not your source files include the headers of
#another library, then it is an INTERFACE dependency.
target_link_libraries(GiSmoPlugin LINK_PUBLIC gismo_static)

#get_property(gsParaviewPlugin_location TARGET gsParaviewPlugin PROPERTY LOCATION)

#EXECUTE_PROCESS(COMMAND test -w ${ParaView_DIR}/lib RESULT_VARIABLE canwrite)
#if(${canwrite})
#  ADD_CUSTOM_COMMAND(
#    TARGET GiSmoPlugin
#    POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:GiSmoPlugin> ${ParaView_DIR}/lib
#    COMMAND ${CMAKE_COMMAND} -E echo 'The plugin has been copied to ${ParaView_DIR}/lib.'
#    )
#else()
#  message(WARNING "No write permissions on plugins directory ${ParaView_DIR}")
#endif()
