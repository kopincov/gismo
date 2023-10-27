### CMakeLists.txt ---
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## CFX5 extension
project(gsCFX5Extension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

# Look for CFX5 library
find_package(CFX5 QUIET REQUIRED)

include_directories(${CFX5_INCLUDES})

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
set(gismo_dev_LINKER ${gismo_dev_LINKER} ${CFX5_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsCFX5/
        FILES_MATCHING PATTERN "*.h")
