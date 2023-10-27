######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

project(unittest_dev)

# for CMake >=3.10
#set(CMAKE_DIRECTORY_LABELS unittest_dev)

include(gsConfig)

include_directories(${GSUNITTEST_INCLUDE_DIR})

include(CTest)
enable_testing()

# Collect source file names
aux_gs_cpp_directory(${PROJECT_SOURCE_DIR} unittest_dev_SRCS)

get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} PATH)
add_definitions(-DGISMO_DDATA_DIR="${PARENT_DIR}/filedata/")

add_executable(unittest_dev gismo_unittest.h ${unittest_dev_SRCS} main.cpp)
target_link_libraries(unittest_dev gismo gismo_dev UnitTestPP)
set_property(TARGET unittest_dev PROPERTY LABELS unittest_dev)
set_target_properties(${PROJECT_NAME} PROPERTIES LABELS unittest_dev FOLDER "unittests")

foreach(file ${unittest_dev_SRCS})
    get_filename_component(testname ${file} NAME_WE) # name without extension
    add_test(NAME unit_${testname} COMMAND $<TARGET_FILE:unittest_dev> ${file})
    set_property(TEST unit_${testname} PROPERTY LABELS unittest_dev)
endforeach(file ${FILES})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
