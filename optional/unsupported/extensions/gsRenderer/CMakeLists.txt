### CMakeLists.txt ---
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## Renderer extension
project(gsRendererExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

# OpenGL and GLU library
find_package(OpenGL QUIET)
if (NOT OPENGL_GLU_FOUND)
  add_definitions("-DGLEW_NO_GLU")
endif()

# GLEW
message("STATUS CMake Will attempt to download GLEW sources.")
include(ExternalProject)
ExternalProject_Add(GLEW
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/GLEW
          BUILD_IN_SOURCE 1
          URL https://sourceforge.net/projects/glew/files/glew/1.13.0/glew-1.13.0.tgz
          URL_MD5 7cbada3166d2aadfc4169c4283701066
          PATCH_COMMAND     
                       patch -p0 -N -d ${CMAKE_CURRENT_BINARY_DIR}/GLEW/build/cmake/ < ${CMAKE_CURRENT_SOURCE_DIR}/GLEW.patch
          CMAKE_ARGS
                       -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                       -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                       -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
          CONFIGURE_COMMAND 
                       "${CMAKE_COMMAND}"
                       "${CMAKE_ARGS}"
                       "-G${CMAKE_GENERATOR}"
                       "${CMAKE_CURRENT_BINARY_DIR}/GLEW/build/cmake"
          #BUILD_COMMAND ""
          UPDATE_COMMAND ""
          #INSTALL_COMMAND ""
)

ExternalProject_Get_Property(GLEW install_dir)
set(GLEW_LIBRARY ${install_dir}/lib/libglew.so)
set(GLEW_INCLUDE_DIR ${install_dir}/include CACHE INTERNAL "")
include_directories(${GLEW_INCLUDE_DIR})

# GLM
message("STATUS CMake Will attempt to download GLM sources.")
add_definitions("-DGLM_ENABLE_EXPERIMENTAL")
include(ExternalProject)
ExternalProject_Add(GLM
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/GLM
          URL https://github.com/g-truc/glm/archive/0.9.9.0.tar.gz
          URL_MD5 3e8d6608d43647e1646f94bc5a3edfbe
          CMAKE_ARGS
                       -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                       -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
          #CONFIGURE_COMMAND ""
          #BUILD_COMMAND ""
          UPDATE_COMMAND ""
          #INSTALL_COMMAND ""
)

ExternalProject_Get_Property(GLM install_dir)
set(Glm_DIR ${install_dir}/include CACHE INTERNAL "")
include_directories(${Glm_DIR})

message("STATUS CMake Will attempt to download Renderer sources.")
include(ExternalProject)
ExternalProject_Add(Renderer
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Renderer
          SVN_REPOSITORY https://ricamsvn.ricam.oeaw.ac.at/gismo/renderer
          CMAKE_ARGS
                       -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                       -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                       -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
                       -DGLEW_INCLUDE_DIR=${GLEW_INCLUDE_DIR}
                       -DGLEW_LIBRARY=${GLEW_LIBRARY}
                       -DGlm_DIR=${Glm_DIR}
                       -Dgismo_DIR=${CMAKE_BINARY_DIR}
          #CONFIGURE_COMMAND ""
          #BUILD_COMMAND ""
          UPDATE_COMMAND ""
          INSTALL_COMMAND ""
          )

ExternalProject_Get_Property(Renderer source_dir)
set(RENDERER_INCLUDES "${source_dir}/src;${Glm_DIR};${GLEW_INCLUDE_DIR}" CACHE INTERNAL "")
include_directories(${RENDERER_INCLUDES})

# Qt includes
find_package(Qt5 COMPONENTS Widgets OpenGL QUIET)
if (Qt5_FOUND)
  message(STATUS "Found Qt5")
  include_directories(${Qt5Widgets_INCLUDE_DIRS}
                      ${Qt5OpenGL_INCLUDE_DIRS})
#  qt5_use_modules(Renderer Core Gui Widgets OpenGL)
endif()

set (QT_USE_QTXML 1)
set (QT_USE_QTOPENGL 1)
add_definitions (${QT_DEFINITIONS})
add_definitions (-DQT_PLUGIN)
add_definitions (-DQT_SHARED)
add_definitions (-DQT_NO_DEBUG)

add_dependencies(Renderer gismo)
add_dependencies(Renderer GLM)
add_dependencies(Renderer GLEW)

add_library(${PROJECT_NAME} OBJECT
  ${${PROJECT_NAME}_SOURCES}
  ${${PROJECT_NAME}_HEADERS}
  ${${PROJECT_NAME}_HPPFILES}
  )

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
  FOLDER "G+Smo extensions"
  )

set(gismo_dev_EXTENSIONS ${gismo_dev_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "${PROJECT_NAME} extensions to be included")

add_dependencies(${PROJECT_NAME} Renderer)

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsRenderer/
        FILES_MATCHING PATTERN "*.h")

#ExternalProject_Get_Property(Renderer source_dir)
#set(RENDERER_INCLUDES "${source_dir}/src;${Glm_DIR};${GLEW_INCLUDE_DIR}" CACHE INTERNAL "")
#include_directories(${RENDERER_INCLUDES})

#ExternalProject_Get_Property(Renderer install_dir)
#set(RENDERER_LIBRARIES ${install_dir}/src/Renderer-build/lib/librenderer.so CACHE INTERNAL "")
