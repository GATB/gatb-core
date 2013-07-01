project(org.inria.gatb.tools)

cmake_minimum_required (VERSION 2.6)

################################################################################
# Define cmake modules directory
################################################################################
set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


################################################################################
# Include external cmake files
################################################################################
include (DefineInteger)
include (CppUnit)

DefineInteger (k)

################################################################################
# OPENMP MANAGEMENT
################################################################################
find_package(OpenMP)

if (OPENMP_FOUND AND NOT (omp EQUAL 0))
    set (COMPILE_DEFINITIONS    "${COMPILE_DEFINITIONS}    ${OpenMP_C_FLAGS} -DOMP")
endif()

################################################################################
# GENERAL DEFINITIONS
################################################################################
if (debug)
    set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS} -g -p -pg")
    message("-- COMPILATION IN DEBUG MODE")
else()
    set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS} -O3")
endif()

set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS}  -DINTEGER_KIND=${INTEGER_KIND}  -DKMER_PRECISION=${KMER_PRECISION}")

if (INT128_FOUND)
    set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS}  -DINT128_FOUND")
endif()

message("-- Options: ${COMPILE_DEFINITIONS}")

add_definitions (${COMPILE_DEFINITIONS})

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
set (gatb-core-flags ${COMPILE_DEFINITIONS})

# We define the include directories used for linking binary based on gatb core
set (gatb-core-includes   ${CMAKE_CURRENT_SOURCE_DIR}/src  ${CMAKE_CURRENT_SOURCE_DIR}/src/gatb/tools)

# We define the libraries used for linking binary based on gatb core
set (gatb-core-libraries   GatbTools-static  dl  pthread  z  ${OpenMP_C_FLAGS})

# NOTE... we have to duplicate the variables for the other scopes (in particular for sub directories)
set (gatb-core-flags     ${gatb-core-flags}     PARENT_SCOPE)
set (gatb-core-includes  ${gatb-core-includes}  PARENT_SCOPE)
set (gatb-core-libraries ${gatb-core-libraries} PARENT_SCOPE)

################################################################################
#  UNIT TESTS GENERATION 
################################################################################
IF (DEFINED CPPUNIT_FOUND)
    ADD_SUBDIRECTORY(test )
ENDIF()

################################################################################
#  DOCUMENTATION GENERATION 
################################################################################
ADD_SUBDIRECTORY(doc EXCLUDE_FROM_ALL)

################################################################################
#  EXAMPLES GENERATION 
################################################################################
ADD_SUBDIRECTORY(examples EXCLUDE_FROM_ALL)

################################################################################
#  TOOLS GENERATION 
################################################################################
ADD_SUBDIRECTORY(tools EXCLUDE_FROM_ALL)