######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2012-2018 - RICAM-Linz.
######################################################################

#if("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
#   message(FATAL_ERROR "The devel folder cannot be configured alone, configure /trunk instead. Remove garbage CMake files and try again.")
#endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

cmake_minimum_required(VERSION 2.8.8)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 NEW)
cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0048 NEW)
endif()

if(POLICY CMP0054)# CMake 3.1
cmake_policy(SET CMP0054 NEW)
endif()

if(POLICY CMP0074)# CMake 3.12
cmake_policy(SET CMP0074 NEW)
endif()

if(POLICY CMP0063)# CMake 3.3
 cmake_policy(SET CMP0063 NEW)
endif()

set(gismo_dev_VERSION_MAJOR "${gismo_VERSION_MAJOR}")
set(gismo_dev_VERSION_MINOR "${gismo_VERSION_MINOR}")
set(gismo_dev_VERSION_PATCH "${gismo_VERSION_PATCH}")
set(gismo_dev_VERSION
  "${gismo_dev_VERSION_MAJOR}.${gismo_dev_VERSION_MINOR}.${gismo_dev_VERSION_PATCH}" CACHE INTERNAL "gismo_dev version number")

if(POLICY CMP0048)# CMake 3.0
    project(gismo_dev LANGUAGES C CXX  VERSION ${gismo_VERSION})
else()
    project(gismo_dev C CXX)
endif()

if (GISMO_BUILD_PCH)
    enable_language(CXXPCH)
endif()

if(NOT GISMO_INSOURCE_BUILD AND "${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source G+Smo builds are not permitted. \n How to fix this:\n 1. Type \nmake cleancmake\nto remove garbage CMake files. \n\n 2. Please make a separate folder for building, otherwise typing \nmake\nthat will create a ./build folder and will compile in that folder. ")
endif()

## #################################################################
## Configuration and Compilation
## #################################################################

message(STATUS "Configuring Developer's branch")
include( gsConfig )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include( gsDevOptions )

## #################################################################
## Preprocessor definitions
## #################################################################

# Shared pointer
find_package (TR1 QUIET)

# Data directory for devel
set(GISMO_DEV_DATA_DIR "${PROJECT_SOURCE_DIR}/filedata/")

# Configure gismo_dev DLL exporting
configure_file ("${PROJECT_SOURCE_DIR}/src/gsCore/gsDevConfig.h.in"
                "${PROJECT_BINARY_DIR}/gsCore/gsDevConfig.h" )
# file(COPY "${PROJECT_SOURCE_DIR}/CTestConfig.cmake"
#      DESTINATION "${PROJECT_BINARY_DIR}")

## #################################################################
## Collect headers and sources
## #################################################################

# Clean previously set variables in the cache
unset(${PROJECT_NAME}_HEADERS    CACHE)
unset(${PROJECT_NAME}_SOURCES    CACHE)
unset(${PROJECT_NAME}_EXTENSIONS CACHE)
unset(${PROJECT_NAME}_MODULES    CACHE)
unset(${PROJECT_NAME}_examples   CACHE)
unset(${PROJECT_NAME}_optionals  CACHE)
unset(${PROJECT_NAME}_LINKER     CACHE)
unset(gismo_TMPL_HEADERS CACHE)
set(${PROJECT_NAME}_optionals ""
    CACHE INTERNAL "${PROJECT_NAME} list of optionals")

# Add main sources directory
set (GISMO_DEV_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src;${PROJECT_SOURCE_DIR}/external;${PROJECT_SOURCE_DIR}/extensions;${PROJECT_BINARY_DIR}" CACHE INTERNAL "${PROJECT_NAME} include dirs")

include_directories(${GISMO_INCLUDE_DIRS})
include_directories(${GISMO_DEV_INCLUDE_DIRS})
include_directories(${gismo_BINARY_DIR}) # For trunk !
add_subdirectory(src)

# Add externals directory
#add_subdirectory(external)

## #################################################################
## Extensions
## #################################################################

if (GISMO_WITH_CAPI)
   add_subdirectory(${PROJECT_SOURCE_DIR}/extensions/gsCApi)
endif()

if(GISMO_WITH_CFX5)
  add_subdirectory(extensions/gsCFX5)
  set (GISMO_DEV_INCLUDE_DIRS ${GISMO_DEV_INCLUDE_DIRS} ${CFX5_INCLUDES}
    CACHE INTERNAL "${PROJECT_NAME} include directories")
  set(gismo_dev_LINKER ${gismo_dev_LINKER} ${CFX5_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  include_directories(${CFX5_INCLUDES})
endif(GISMO_WITH_CFX5)

if (GISMO_WITH_RENDERER)
   add_subdirectory(${PROJECT_SOURCE_DIR}/extensions/gsRenderer)
   set (GISMO_DEV_INCLUDE_DIRS ${GISMO_DEV_INCLUDE_DIRS} ${RENDERER_INCLUDES}
     CACHE INTERNAL "${PROJECT_NAME} include directories")
  set(gismo_dev_LINKER ${gismo_dev_LINKER} ${RENDERER_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
   include_directories(${RENDERER_INCLUDES})
endif()

if (GISMO_WITH_VTK)
  add_subdirectory(extensions/gsVtk)
endif(GISMO_WITH_VTK)

## #################################################################
## Find stable branch
## #################################################################

#find_package(gismo REQUIRED CMAKE_FIND_ROOT_PATH_BOTH)
#link_directories   ( ${GISMO_LIBRARY_DIR})
#include( ${GISMO_USE_FILE})

## #################################################################
## Build library instantized for GISMO_COEFF_TYPE number type
## #################################################################

include(gsLibrary)

add_custom_target(unsupported) #alias
add_dependencies(unsupported ${PROJECT_NAME})

if(GISMO_BUILD_LIB)
  target_link_libraries(${PROJECT_NAME} gismo)
endif()

#message ("*** Headers:  ${${PROJECT_NAME}_HEADERS}")
#message ("*** Sources:  ${${PROJECT_NAME}_SOURCES}")
#message ("*** Template headers:  ${${PROJECT_NAME}_TMPL_HEADERS}")


## #################################################################
## Plugins and optionals
## #################################################################

add_custom_target(plugins_dev)

add_subdirectory(optional)

if(GISMO_BUILD_PVIEW)
  add_subdirectory(plugins/gsParaview)
endif(GISMO_BUILD_PVIEW)

## #################################################################
## Configure / install
## #################################################################

#include( gsInstall-dev )
include(${PROJECT_SOURCE_DIR}/cmake/gsInstall-dev.cmake)

## #################################################################
## Executables
## #################################################################

if(GISMO_BUILD_UNITTESTS)
  add_subdirectory(unittests)
else()
  add_custom_target(unittest_dev)
endif(GISMO_BUILD_UNITTESTS)

## #################################################################
## Examples
## #################################################################

if(GISMO_BUILD_EXAMPLES)
  add_subdirectory(examples)
else()
  add_subdirectory(examples EXCLUDE_FROM_ALL)
endif(GISMO_BUILD_EXAMPLES)

## #################################################################
## Documentation
## #################################################################

add_subdirectory(doc)
