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
    #set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS} -g -p -pg")
    set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS} -g -p")
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
#  DIRECTORIES MANAGEMENT 
################################################################################

set (LIBRARY_OUTPUT_PATH       ${PROJECT_BINARY_DIR}/lib/${CMAKE_BUILD_TYPE})
set (EXECUTABLE_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE})

################################################################################
#  EXPORT FOR THIRD PARTIES 
################################################################################

# We define the compilation flags used for compiling binary based on gatb core
set (gatb-core-flags ${COMPILE_DEFINITIONS})

# We define the include directories used for linking binary based on gatb core
set (gatb-core-includes   
    ${PROJECT_SOURCE_DIR}/src  
    ${PROJECT_SOURCE_DIR}/src/gatb/tools  
    ${PROJECT_SOURCE_DIR}/thirdparty/hdf5/src
    ${PROJECT_BINARY_DIR}/thirdparty/hdf5  
)

# We define the libraries used for linking binary based on gatb core
set (gatb-core-libraries   gatbcore-static  dl  pthread  z hdf5  ${OpenMP_C_FLAGS} )

# NOTE... we have to duplicate the variables for the other scopes (in particular for sub directories)
set (gatb-core-flags     ${gatb-core-flags}     PARENT_SCOPE)
set (gatb-core-includes  ${gatb-core-includes}  PARENT_SCOPE)
set (gatb-core-libraries ${gatb-core-libraries} PARENT_SCOPE)

################################################################################
#  LIBRARY GENERATION 
################################################################################
ADD_SUBDIRECTORY(src)

################################################################################
#  UNIT TESTS GENERATION 
################################################################################
IF (DEFINED CPPUNIT_FOUND)
    ADD_SUBDIRECTORY (test EXCLUDE_FROM_ALL)
ENDIF()

################################################################################
#  TOOLS GENERATION 
################################################################################
ADD_SUBDIRECTORY(tools)

################################################################################
#  THIRD PARTY GENERATION (
################################################################################
ADD_SUBDIRECTORY(thirdparty)

################################################################################
#  DEPENDENCIES 
################################################################################
ADD_DEPENDENCIES (gatbcore-static hdf5)

################################################################################
#  DOCUMENTATION GENERATION 
################################################################################
ADD_SUBDIRECTORY(doc EXCLUDE_FROM_ALL)

################################################################################
#  EXAMPLES GENERATION 
################################################################################
ADD_SUBDIRECTORY(examples EXCLUDE_FROM_ALL)

