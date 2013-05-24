project(org.inria.gatb.tools)

cmake_minimum_required (VERSION 2.6)


################################################################################
# GENERAL DEFINITIONS
################################################################################
add_definitions (-O3 -funroll-loops  -fomit-frame-pointer  -msse3  -DINTEGER_KIND=1  -DKMER_PRECISION=3 -fopenmp)
#add_definitions (-g -p -pg  -DKMER_PRECISION=2)


################################################################################
# OPERATING SYSTEM SPECIFIC
################################################################################
IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions(-D__WINDOWS__)
    # IMPORTANT ! we have to add this flag for getting correct methods names in the library    
    # otherwise the Java JVM won't find the native JNI methods.
    set_target_properties (PlastLibrary PROPERTIES LINK_FLAGS -Wl,--kill-at)
    set_target_properties (PlastLibrary PROPERTIES PREFIX "")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")


################################################################################
#  LIBRARY GENERATION 
################################################################################

set (LIBRARY_OUTPUT_PATH    lib/${CMAKE_BUILD_TYPE})

file (GLOB_RECURSE  ComponentFiles  src/*)

include_directories (src  src/gatb/tools)

add_library    (GatbTools-static  STATIC  ${ComponentFiles} )
#add_library    (GatbTools-dynamic SHARED  ${ComponentFiles} )

set_target_properties (GatbTools-static   PROPERTIES OUTPUT_NAME  GatbTools  clean_direct_output 1)
#set_target_properties (GatbTools-dynamic  PROPERTIES OUTPUT_NAME  GatbTools  clean_direct_output 1)


################################################################################
#  UNIT TESTS GENERATION 
################################################################################
ADD_SUBDIRECTORY(test)

################################################################################
#  DOCUMENTATION GENERATION 
################################################################################
ADD_SUBDIRECTORY(doc)

################################################################################
#  EXAMPLES GENERATION 
################################################################################
ADD_SUBDIRECTORY(examples)

################################################################################
#  TOOLS GENERATION 
################################################################################
ADD_SUBDIRECTORY(tools)
