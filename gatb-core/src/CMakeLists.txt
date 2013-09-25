################################################################################
#  LIBRARY GENERATION 
################################################################################

file (GLOB_RECURSE  ComponentFiles  *)

# WARNING !!! For the moment, we need to remove some warnings (on Macos) due to use of offsetof macro on non Plain Old Data
set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS}  -Wno-invalid-offsetof")

include_directories (${gatb-core-includes})

add_library    (gatbcore-static  STATIC  ${ComponentFiles} )
#add_library   (gatbcore-dynamic SHARED  ${ComponentFiles} )

set_target_properties  (gatbcore-static   PROPERTIES OUTPUT_NAME  gatbcore  clean_direct_output 1)
#set_target_properties (gatbcore-dynamic  PROPERTIES OUTPUT_NAME  gatbcore  clean_direct_output 1)
