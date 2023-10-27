## #################################################################
## Installation
## #################################################################

#message ("  CMAKE_INSTALL_PREFIX    ${CMAKE_INSTALL_PREFIX}")

#set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR     lib     CACHE STRING "Installation directory for libraries")
set(INSTALL_BIN_DIR     bin     CACHE STRING "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE STRING "Installation directory for header files")

# Set CMake installation directory
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR ${INSTALL_LIB_DIR}/cmake)
else()
  #set(DEF_INSTALL_CMAKE_DIR ${INSTALL_LIB_DIR}/cmake/gismo)
   set(DEF_INSTALL_CMAKE_DIR ${INSTALL_LIB_DIR})
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# Add all targets to the build-tree export set
if(GISMO_BUILD_LIB)
export(TARGETS ${PROJECT_NAME}
  FILE "${PROJECT_BINARY_DIR}/gismoTargets.cmake" APPEND)
endif()

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE gismo)

# Create the gismo_devConfig.cmake

# ... for the build tree
set(CONF_INCLUDE_DIRS "${GISMO_DEV_INCLUDE_DIRS}"
                      "${PROJECT_SOURCE_DIR}/src"
                      "${PROJECT_SOURCE_DIR}/external"
                      "${PROJECT_BINARY_DIR}"
                      #"${PROJECT_SOURCE_DIR}/extensions"
		      )
set(CONF_LIB_DIRS     "${CMAKE_BINARY_DIR}/lib")
configure_file(${PROJECT_SOURCE_DIR}/cmake/gismo_devConfig.cmake.in
  "${CMAKE_BINARY_DIR}/gismo_devConfig.cmake" @ONLY)

# ... for the install tree
set(CONF_INCLUDE_DIRS "${INSTALL_INCLUDE_DIR}/${PROJECT_NAME}")
set(CONF_LIB_DIRS     "${INSTALL_LIB_DIR}")
configure_file(${PROJECT_SOURCE_DIR}/cmake/gismo_devConfig.cmake.in
  "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/gismo_devConfig.cmake" @ONLY)

if(GISMO_BUILD_LIB)

# For gsExport.h (upd: now shared with stable)
#install(FILES ${PROJECT_BINARY_DIR}/gsCore/gsExport.h DESTINATION include/${PROJECT_NAME}/gsCore/)

# Install cmake files
install(FILES
  "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/gismo_devConfig.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

else(GISMO_BUILD_LIB)
   message ("Configure with -DGISMO_BUILD_LIB=ON to compile the library")
endif(GISMO_BUILD_LIB)
