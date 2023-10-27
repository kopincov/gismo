######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

project(example_dev)

# Add a grouping target that builds all tests this is equivalent to
# running "make" inside the tests build folder
add_custom_target(${PROJECT_NAME})
set_target_properties(${PROJECT_NAME} PROPERTIES LABELS "${PROJECT_NAME}" FOLDER "${PROJECT_NAME}-gismo")

include(CTest)
enable_testing()

# Collect source file names

aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR} FILES)

if(NOT GISMO_BUILD_CPPLOT)
  list(REMOVE_ITEM FILES ${CMAKE_CURRENT_SOURCE_DIR}/gsPlot_test.cpp)
endif(NOT GISMO_BUILD_CPPLOT)

foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
	set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}-gismo")
    add_dependencies(${PROJECT_NAME} ${tarname})
    if (GISMO_BUILD_LIB)
      target_link_libraries(${tarname} gismo_dev)
    else()
      target_link_libraries(${tarname} gismo_dev_static)
    endif()
endforeach(file ${FILES})

if(NOT WIN32)
target_link_libraries(gsJITCompiler_test dl)
endif()

list(LENGTH FILES len)
message(STATUS "Number of tests to compile: ${len}")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
