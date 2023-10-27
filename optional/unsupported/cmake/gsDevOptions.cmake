######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2016 RICAM-Linz.
######################################################################

## #################################################################
## Options list: Dev options
## #################################################################

option(GISMO_WITH_CFX5           "With CApi"                     false  )
if (${GISMO_WITH_CAPI})
message ("  GISMO_WITH_CAPI         ${GISMO_WITH_CAPI}")
endif()

option(GISMO_WITH_CAPI           "With CFX5"                     false  )
if (${GISMO_WITH_CFX5})
message ("  GISMO_WITH_CFX5         ${GISMO_WITH_CFX5}")
endif()

option(GISMO_WITH_RENDERER       "With Renderer"                 false  )
if (${GISMO_WITH_RENDERER})
message ("  GISMO_WITH_RENDERER     ${GISMO_WITH_RENDERER}")
endif()
