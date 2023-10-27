######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

cmake_minimum_required(VERSION 2.8.8)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 NEW)
cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0048 NEW)
endif()

if(POLICY CMP0054)# CMake 3.1
cmake_policy(SET CMP0054 NEW)
endif()

if(POLICY CMP0063)# CMake 3.3
cmake_policy(SET CMP0063 NEW)
endif()

project(motor)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_custom_target(${PROJECT_NAME})
set_property(TARGET ${PROJECT_NAME} PROPERTY LABELS "${PROJECT_NAME}")

include(gsConfig)

unset(${PROJECT_NAME}_HEADERS CACHE)
#message(" {GISMO_DEV_INCLUDE_DIRS} ${GISMO_DEV_INCLUDE_DIRS} ")

include_directories(${GISMO_INCLUDE_DIRS})
include_directories(${GISMO_DEV_INCLUDE_DIRS})
include_directories(${gismo_BINARY_DIR}) # For trunk !
include_directories(src)

add_definitions(-DMOTOR_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/filedata/")

add_subdirectory(jku)
add_subdirectory(tud)
add_subdirectory(uwb-pilsen)

include(CTest)
enable_testing()

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
set(LIBRARY_OUTPUT_PATH    ${CMAKE_BINARY_DIR}/lib/)
