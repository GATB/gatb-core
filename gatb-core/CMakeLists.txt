project(gatb-core)

# We set the required version
cmake_minimum_required (VERSION 2.6)

################################################################################
# The version number.
################################################################################
SET (gatb-core_VERSION_MAJOR 1)
SET (gatb-core_VERSION_MINOR 0)
SET (gatb-core_VERSION_PATCH 0)

IF (DEFINED MAJOR)
    SET (gatb-core_VERSION_MAJOR ${MAJOR})
ENDIF()
IF (DEFINED MINOR)
    SET (gatb-core_VERSION_MINOR ${MINOR})
ENDIF()
IF (DEFINED PATCH)
    SET (gatb-core_VERSION_PATCH ${PATCH})
ENDIF()

set (gatb-core-version ${gatb-core_VERSION_MAJOR}.${gatb-core_VERSION_MINOR}.${gatb-core_VERSION_PATCH})
set (gatb-core-date "xxxx-xx-xx")

################################################################################
# Boost
################################################################################

# We need boost for its "variant" part. It could be worthy to estimate a custom
# (and possibly optimized) version of our own in order to avoid the boost dependency.
find_package (Boost REQUIRED)

################################################################################
# Define cmake modules directory
################################################################################
set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

################################################################################
# Include external cmake files
################################################################################
include (DefineInteger)
include (CppUnit)
include (GetDate)

# We check whether we have native 128 bits integers
DefineInteger (k)

# We get the current date
GetCurrentDate (gatb-core-date) 

################################################################################
# GENERAL DEFINITIONS
################################################################################
if (debug)
    #set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS} -g -p -pg")
    set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS} -g -p")
    message("-- COMPILATION IN DEBUG MODE")
else()
    set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS} -O3")
endif()

if (INT128_FOUND)
    set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS}  -DINT128_FOUND")
endif()

# We determine what kind of c++ standard we should use according to the compiler
# This should work for gcc and clang
set (CXX_STD_VERSION "c++11")

# In case we use an "old" version of gcc, we can still do something with the c++0x standard
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU"  AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7")
    set (CXX_STD_VERSION "c++0x")
endif()

# We use lambda expressions, so we need to use the correct compiler option
set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS}  -std=${CXX_STD_VERSION}")

message("-- Options: ${LIBRARY_COMPILE_DEFINITIONS}")

################################################################################
#  DIRECTORIES MANAGEMENT 
################################################################################

set (LIBRARY_OUTPUT_PATH       ${PROJECT_BINARY_DIR}/lib/${CMAKE_BUILD_TYPE})
set (EXECUTABLE_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE})

################################################################################
#  EXPORT FOR THIRD PARTIES 
################################################################################

# We define the compilation flags used for compiling binary based on gatb core
set (gatb-core-flags ${LIBRARY_COMPILE_DEFINITIONS})

# We define the include directories used for linking binary based on gatb core
# Note that we need to add boost include dependency (not working otherwise with clang)
set (gatb-core-includes ${PROJECT_BINARY_DIR}/include  ${Boost_INCLUDE_DIRS})

# We define the libraries used for linking binary based on gatb core
set (gatb-core-libraries   gatbcore-static  dl  pthread  z hdf5)

# We define the directory where to find cmake helpers
set (gatb-core-cmake  ${CMAKE_MODULE_PATH})

# NOTE... we have to duplicate the variables for the other scopes (in particular for sub directories)
set (gatb-core-flags     ${gatb-core-flags}     PARENT_SCOPE)
set (gatb-core-includes  ${gatb-core-includes}  PARENT_SCOPE)
set (gatb-core-libraries ${gatb-core-libraries} PARENT_SCOPE)
set (gatb-core-cmake     ${gatb-core-cmake}     PARENT_SCOPE)

################################################################################
#  LIBRARY GENERATION 
################################################################################
ADD_SUBDIRECTORY(src)

################################################################################
#  UNIT TESTS GENERATION 
################################################################################
IF (DEFINED CPPUNIT_FOUND)
    ADD_SUBDIRECTORY (test)
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
# we must be sure that hdf5 is built and installed before building gatb-core
ADD_DEPENDENCIES (gatbcore-static hdf5 hdf5_postbuild)

################################################################################
#  DOCUMENTATION GENERATION 
################################################################################
ADD_SUBDIRECTORY(doc EXCLUDE_FROM_ALL)

################################################################################
#  EXAMPLES GENERATION 
################################################################################
ADD_SUBDIRECTORY(examples EXCLUDE_FROM_ALL)

################################################################################
#  INSTALL 
################################################################################
IF (NOT DEFINED GATB_CORE_INSTALL_EXCLUDE)
    INSTALL (FILES ${PROJECT_SOURCE_DIR}/doc/misc/README.txt  DESTINATION . OPTIONAL)
    INSTALL (DIRECTORY ${PROJECT_SOURCE_DIR}/examples DESTINATION . OPTIONAL  FILES_MATCHING PATTERN "*.cpp")
ENDIF()

################################################################################
#  DELIVERY 
################################################################################
SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "gatb-core project")
SET (CPACK_PACKAGE_VENDOR               "Genscale team (INRIA)")
SET (CPACK_PACKAGE_VERSION_MAJOR        "${gatb-core_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR        "${gatb-core_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH        "${gatb-core_VERSION_PATCH}")
SET (CPACK_PACKAGE_VERSION              "${gatb-core-version}")  
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")
SET (CPACK_SOURCE_IGNORE_FILES          
    "^${PROJECT_SOURCE_DIR}/build/"  
    "^${PROJECT_SOURCE_DIR}/.project"
    "^${PROJECT_SOURCE_DIR}/.gitignore"
    "^${PROJECT_SOURCE_DIR}/doc/design"
    "^${PROJECT_SOURCE_DIR}/DELIVERY.md"
)

# We include the module and get all the delivery targets
include (Delivery)

