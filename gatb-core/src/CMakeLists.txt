################################################################################
#  LIBRARY GENERATION 
################################################################################

file (GLOB_RECURSE  ComponentFiles  *)

include_directories (${gatb-core-includes})

add_library    (gatbcore-static  STATIC  ${ComponentFiles} )
#add_library   (gatbcore-dynamic SHARED  ${ComponentFiles} )

set_target_properties  (gatbcore-static   PROPERTIES OUTPUT_NAME  gatbcore  clean_direct_output 1)
#set_target_properties (gatbcore-dynamic  PROPERTIES OUTPUT_NAME  gatbcore  clean_direct_output 1)
