### CMakeLists.txt ---
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2013 - RICAM-Linz.
######################################################################

project(gsCApi)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)
#message ("*** ${PROJECT_NAME} sources: ${${PROJECT_NAME}_HEADERS}, ${${PROJECT_NAME}_SOURCES}")

# Look for parasolid library
FIND_PACKAGE(Parasolid REQUIRED)
include_directories(${PARASOLID_INCLUDE_DIR})

# Compile statically the extension
add_library(${PROJECT_NAME} OBJECT
${${PROJECT_NAME}_SOURCES}
${${PROJECT_NAME}_HEADERS}
${${PROJECT_NAME}_HPPFILES}
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE C
    FOLDER "G+Smo extensions"
    )

  set(gismo_dev_MODULES ${gismo_dev_MODULES} $<TARGET_OBJECTS:${PROJECT_NAME}>
      CACHE INTERNAL "G+Smo-dev extensions to be included" )

#if(GISMO_BUILD_LIB)
#install(FILES   ${CMAKE_CURRENT_SOURCE_DIR}/gsCApi.h
#        DESTINATION include/gismo/gsCApi)
#endif()
