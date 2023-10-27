 ### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

#get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} PATH)
add_definitions(-DGISMO_DDATA_DIR="${PROJECT_SOURCE_DIR}/filedata/")


## IpOpt option
if(GISMO_WITH_IPOPT)
  add_gismo_executable(${CMAKE_CURRENT_SOURCE_DIR}/gsOptParam.cpp)
  mark_gismo_optional (gsOptParam.cpp)
  target_link_libraries(gsOptParam ${IPOPT_LIBRARY}) # for header-only

  add_gismo_executable(${CMAKE_CURRENT_SOURCE_DIR}/gsOptMasonry.cpp)
  mark_gismo_optional (gsOptMasonry.cpp)
  target_link_libraries(gsOptMasonry ${IPOPT_LIBRARY})

  add_gismo_executable(${CMAKE_CURRENT_SOURCE_DIR}/gsOptProblemExample.cpp)
  mark_gismo_optional (gsOptProblemExample.cpp)
  target_link_libraries(gsOptProblemExample ${IPOPT_LIBRARY})

  #include_directories(${IPOPT_INCLUDE_DIR})

  ## CoDiPack option
  if(GISMO_WITH_CODIPACK)
    include_directories(${CODIPACK_INCLUDE_DIR})
    if(GISMO_WITH_ADDSL)
      include_directories(${ADDSL_INCLUDE_DIR})
    endif()
    add_gismo_executable(${CMAKE_CURRENT_SOURCE_DIR}/gsOptADtest.cpp)
    mark_gismo_optional (gsOptADtest.cpp)
    target_link_libraries(gsOptADtest ${IPOPT_LIBRARY})
  endif(GISMO_WITH_CODIPACK)

endif(GISMO_WITH_IPOPT)

## CoDiPack option
if(GISMO_WITH_CODIPACK)
  include_directories(${CODIPACK_INCLUDE_DIR})
  if(GISMO_WITH_ADDSL)
    include_directories(${ADDSL_INCLUDE_DIR})
  endif()
  add_gismo_executable(${CMAKE_CURRENT_SOURCE_DIR}/gsCoDiPacktest.cpp)
  mark_gismo_optional (gsCoDiPacktest.cpp)
endif(GISMO_WITH_CODIPACK)

## VTK option
if(GISMO_WITH_VTK)

  add_gismo_executable(gsVtkView_test.cpp)
  mark_gismo_optional(gsVtkView_test.cpp)

endif(GISMO_WITH_VTK)


### cppplot option
#if(GISMO_WITH_CPPLOT)
## Link with cpplot for plotting
#target_link_libraries(gsPlot_test 
#cpplot
##boostpart
#)
#endif(GISMO_BUILD_CPPLOT)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
