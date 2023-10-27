######################################################################
## CMakeLists.txt ---
## This file is NOT part of the G+Smo library.
##
## Author: MOTOR
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

include(CTest)
enable_testing()

set(SUBPROJECT_NAME motor-tud)

add_custom_target(${SUBPROJECT_NAME})
set_property(TARGET ${PROJECT_NAME} PROPERTY LABELS "${PROJECT_NAME}")

# Collect source file names
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} FILES)
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} HEADERS)
set(${PROJECT_NAME}_HEADERS ${${PROJECT_NAME}_HEADERS} ${HEADERS}
      CACHE INTERNAL "Motor files" )

foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
#    if (GISMO_BUILD_LIB)
      target_link_libraries(${tarname} gismo gismo_dev)
#    endif()

    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME};${SUBPROJECT_NAME}")
    add_dependencies(${PROJECT_NAME}    ${tarname})
    add_dependencies(${SUBPROJECT_NAME} ${tarname})

endforeach(file ${FILES})

list(LENGTH FILES len)
message(STATUS "${SUBPROJECT_NAME} programs to compile: ${len}")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
