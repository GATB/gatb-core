project(org.inria.gatb.tools)

cmake_minimum_required (VERSION 2.6)

################################################################################
# Define cmake modules directory
################################################################################
set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


################################################################################
# DEFINE INTEGER CLASS TO BE USED
################################################################################
include (DefineInteger)
DefineInteger (k)

################################################################################
# OPENMP MANAGEMENT
################################################################################
find_package(OpenMP)

if (OPENMP_FOUND)
    set (COMPILE_DEFINITIONS    "${COMPILE_DEFINITIONS}    ${OpenMP_C_FLAGS} -DOMP")
endif()

################################################################################
# GENERAL DEFINITIONS
################################################################################
if (debug)
    set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS} -g -p -pg")
    message("-- COMPILATION IN DEBUG MODE")
else()
    set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS} -O3 -funroll-loops  -fomit-frame-pointer")
endif()

set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS}  -DINTEGER_KIND=${INTEGER_KIND}  -DKMER_PRECISION=${KMER_PRECISION}")

if (INT128_FOUND)
    set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS}  -DINT128_FOUND")
endif()

message("-- Options: ${COMPILE_DEFINITIONS}")

add_definitions (${COMPILE_DEFINITIONS})

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
#  EXPORT FOR THIRD PARTIES 
################################################################################

# We define the compilation flags used for compiling binary based on gatb core
set (gatb-core-flags ${COMPILE_DEFINITIONS} PARENT_SCOPE)

# We define the include directories used for linking binary based on gatb core
set (gatb-core-includes   ${CMAKE_CURRENT_SOURCE_DIR}/src  ${CMAKE_CURRENT_SOURCE_DIR}/src/gatb/tools  PARENT_SCOPE)

# We define the libraries used for linking binary based on gatb core
set (gatb-core-libraries   GatbTools-static  dl  pthread  z  ${OpenMP_C_FLAGS}  PARENT_SCOPE)

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