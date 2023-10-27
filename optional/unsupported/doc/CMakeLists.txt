######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2012-2015 - RICAM-Linz.
######################################################################

project(GismoDoc)

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

find_package(Doxygen)

if(DOXYGEN_FOUND)

project(GismoDoc-dev)

if(CMAKE_COMPILER_IS_GNUCXX)
  if(CMAKE_SYSTEM_NAME MATCHES Linux)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -g1")
  endif(CMAKE_SYSTEM_NAME MATCHES Linux)
endif(CMAKE_COMPILER_IS_GNUCXX)

# Set some Doxygen flags
set(GISMO_DOXY_PROJECT_NAME             "G+Smo")
set(GISMO_DOXY_OUTPUT_DIRECTORY_SUFFIX  "")
set(GISMO_DOXY_INPUT                    "\"${gismo_SOURCE_DIR}/src\" \"${gismo_SOURCE_DIR}/doc\" \"${gismo_SOURCE_DIR}/examples\" \"${gismo_dev_SOURCE_DIR}/src\" \"${gismo_dev_SOURCE_DIR}/doc\" ")
# Uncomment the line below to have doxygen process only
# what is inside the doc folder when typing 'make doc-devel'
#set(GISMO_DOXY_INPUT                    "\"${gismo_SOURCE_DIR}/doc\" \"${gismo_dev_SOURCE_DIR}/doc\" ")
set(IMAGE_PATH        "\"${gismo_SOURCE_DIR}/doc/figs\" \"${gismo_dev_SOURCE_DIR}/doc/figs\"")
set(GISMO_DOXY_HTML_COLORSTYLE_HUE      "120") #220
set(GISMO_DOXY_HTML_COLORSTYLE_SAT      "100") #100
set(GISMO_DOXY_HTML_COLORSTYLE_GAMMA    "100")  # 80
set(GISMO_DOXY_TAGFILES                 "")
set(SHOW_INCLUDE_FILES                  "yes")
set(CLASS_FLAGS                         "yes")
set(COLLABORATION_GRAPH                 "yes")
set(SHOW_USED_FILES                     "no")
set(INCLUDED_BY_GRAPH                   "no")
set(GROUP_GRAPHS                        "no")
set(UML_LOOK                            "no")

# get the date in variable TODAY
if(CMAKE_MAJOR_VERSION EQUAL 2 AND
   CMAKE_MINOR_VERSION LESS  9 AND
   CMAKE_PATCH_VERSION LESS 11 )
 if (WIN32)
   execute_process(COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE TODAY)
   string(REGEX REPLACE "(..)/(..)/..(..).*" "\\1/\\2/\\3" TODAY ${TODAY})
 elseif(UNIX)
   execute_process(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE TODAY)
   string(REGEX REPLACE "(..)/(..)/..(..).*" "\\1/\\2/\\3" TODAY ${TODAY})
 else (WIN32)
   message(SEND_ERROR "date not implemented")
   set(${RESULT} 000000)
 endif (WIN32)
else()
  string(TIMESTAMP TODAY "%d/%m/%Y")
endif()

#Set default MathJax url
if (NOT MATHJAX_RELPATH)
  set(MATHJAX_RELPATH http://people.ricam.oeaw.ac.at/gismo/mj)
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/gismodoxy_layout.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/gismodoxy_layout.xml COPYONLY)

#add_definitions("-DGISMO_MAKING_DOCS")

set(snippets_targets "")
#add_subdirectory(snippets)

add_custom_target(doc-devel #ALL
  COMMAND doxygen
  WORKING_DIRECTORY ${gismo_dev_BINARY_DIR}/doc
  COMMENT "Generating doc: ${gismo_dev_BINARY_DIR}/doc/html/index.html" VERBATIM
)

#add_custom_command(
#      TARGET doc-devel
#      POST_BUILD
#      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${gismo_dev_BINARY_DIR}/doc/html/usergroup0.html ${gismo_dev_BINARY_DIR}/doc/html/usergroup1.html
#)


# need snippets
#add_dependencies(doc-devel doc)

else (DOXYGEN_FOUND)

add_custom_target(doc-devel
    ${CMAKE_COMMAND} -E cmake_echo_color --cyan "Doxygen was not found! Please install doxygen and run cmake again to build the doc.")

endif(DOXYGEN_FOUND)
