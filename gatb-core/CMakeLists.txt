project(gatb-core)

# We set the required version
cmake_minimum_required (VERSION 2.6)

################################################################################
# The version number.
################################################################################
# The default version number is the latest official build
SET (gatb-core_VERSION_MAJOR 1)
SET (gatb-core_VERSION_MINOR 2)
SET (gatb-core_VERSION_PATCH 1)

# But, it is possible to define another release number during a local build
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

# However, continuous integration has priority over local compilation
IF (DEFINED JENKINS_TAG)
    SET (gatb-core-version ${JENKINS_TAG})
ENDIF()

# User login to upload doc-api; Jenkins has priority over local user
IF (DEFINED JENKINS_GFORGE_USER)
    SET (gatb-doc-user-login ${JENKINS_GFORGE_USER})
ELSE()
    SET (gatb-doc-user-login $ENV{USER})
ENDIF()

set (gatb-core-date "xxxx-xx-xx")

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
# COMPILER DEFINITIONS
################################################################################
# We try to find the best compiler config, in particular for C++ extension usage
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

    IF(NOT CMAKE_CXX_COMPILER_VERSION) #for old cmakes
        include(CheckCompiler)
    ENDIF()
    if (CMAKE_CXX_COMPILER_VERSION  VERSION_LESS 4.5) 
    else()
        set (use_new_cxx 1)
    endif()

    if (CMAKE_CXX_COMPILER_VERSION  VERSION_LESS 4.7) 
        set (use_mphf 0)
    else()
        set (use_mphf 1)
    endif()

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

   IF(CMAKE_SYSTEM_NAME MATCHES "(Darwin)") # different clang versions number between linux and mac
        if (CMAKE_CXX_COMPILER_VERSION  VERSION_LESS 4.3) 
            set (use_mphf 0)
        else()
            set (use_mphf 1)
            set (use_new_cxx 1)
        endif()

    else()
        if (CMAKE_CXX_COMPILER_VERSION  VERSION_LESS 3.2) 
            set (use_mphf 0)
        else()
            set (use_new_cxx 1)
            set (use_mphf 1)
        endif()

    endif()
else()
        message ("-------------------------------------------------------------------------------------")
        message ("-- WARNING !!! YOU USE AN UNKNOWN COMPILER...")
        message ("-------------------------------------------------------------------------------------")
endif()

    
if (use_new_cxx)
    message ("-------------------------------------------------------------------------------------")
    message ("-- C++ EXTENSIONS WILL BE USED...")
    message ("-------------------------------------------------------------------------------------")
else()
    message ("-------------------------------------------------------------------------------------")
    message ("-- COMPILER VERSION IS NOT ENOUGH TO HAVE C++ EXTENSIONS... SOME EXAMPLES CAN'T COMPILE...")
    message ("-------------------------------------------------------------------------------------")
endif()

################################################################################
# GENERAL DEFINITIONS
################################################################################

set (LIB_COMPILE_WARNINGS "-Wall -Wno-unused-function -Wno-unused-variable -Wno-parentheses -Wno-format -Wno-strict-aliasing -Wno-array-bounds -Wno-unknown-pragmas") # last one is exclusively for BooPHF

if (debug)
    #set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS} -g -p -pg")
    set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS} -g -p ${LIB_COMPILE_WARNINGS}")
    set (CMAKE_BUILD_TYPE Debug) # else CMake adds DNDEBUG
    message("-- COMPILATION IN DEBUG MODE")
else()
    set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS} -O3 -DNDEBUG ${LIB_COMPILE_WARNINGS}")
endif()

if (INT128_FOUND)
     set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS}  -DINT128_FOUND")
endif()

if (NONCANONICAL)
    MESSAGE("--- Compiling with non-canonical (direct strand only) k-mer counting")
    set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS}  -DNONCANONICAL=1")
endif()



# We determine what kind of c++ standard we should use according to the compiler
# According to the CXX version, we have to activate some flags
if (use_new_cxx)
    # We set a flag for lambda expressions usage
    set (CXX_STD_VERSION "-std=c++0x")
    # We set a flag for lambda expressions
    set (WITH_LAMBDA_EXPRESSION "1")
    # We add compilation flags to the global flags
    set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS}  -DWITH_LAMBDA_EXPRESSION ${CXX_STD_VERSION}")
endif()

# detect SSE 
# from https://github.com/rurban/smhasher/blob/master/CMakeLists.txt
# i do not see much performance gain for now, but let's keep that code here, might be useful later.
# list of performance gain observed:
# popcount in Graph::countNeighbors
IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
    EXEC_PROGRAM(cat ARGS "/proc/cpuinfo" OUTPUT_VARIABLE CPUINFO)
    STRING(REGEX REPLACE "^.*(sse2).*$" "\\1" SSE_THERE ${CPUINFO})
    STRING(COMPARE EQUAL "sse2" "${SSE_THERE}" SSE2_TRUE)
    STRING(REGEX REPLACE "^.*(sse4_2).*$" "\\1" SSE_THERE ${CPUINFO})
    STRING(COMPARE EQUAL "sse4_2" "${SSE_THERE}" SSE42_TRUE)
ELSEIF(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    EXEC_PROGRAM("/usr/sbin/sysctl -n machdep.cpu.features" OUTPUT_VARIABLE
        CPUINFO)
    STRING(REGEX REPLACE "^.*[^S](SSE2).*$" "\\1" SSE_THERE ${CPUINFO})
    STRING(COMPARE EQUAL "SSE2" "${SSE_THERE}" SSE2_TRUE)
    STRING(REGEX REPLACE "^.*(SSE4.2).*$" "\\1" SSE_THERE ${CPUINFO})
    STRING(COMPARE EQUAL "SSE4.2" "${SSE_THERE}" SSE42_TRUE)
ENDIF(CMAKE_SYSTEM_NAME MATCHES "Linux")
IF (SSE2_TRUE)
    set(SSE2_FOUND true CACHE BOOL "SSE2 available")
ELSE (SSE2_TRUE)
    set(SSE2_FOUND false CACHE BOOL "SSE2 not available")
ENDIF (SSE2_TRUE)
IF (SSE42_TRUE)
    set(SSE4_2_FOUND true CACHE BOOL "SSE4.2 available")
ELSE (SSE42_TRUE)
    set(SSE4_2_FOUND false CACHE BOOL "SSE4.2 not available")
ENDIF (SSE42_TRUE)
# I could use LIBRARY_COMPILE_DEFINITIONS, but it's actually passed to "add_definitions", which isn't made for passing compilation flags, only -D ones.
IF(SSE4_2_FOUND AND (NOT NO_SSE))
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -msse2 -msse4.2 -mpopcnt")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2 -msse4.2 -mpopcnt")
    message ("-- SSE 4.2 detected")
ENDIF()

if (use_new_cxx)
    set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS} -DUSE_NEW_CXX ")
endif()
 
if (use_mphf)
    set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS}    -DWITH_MPHF ")
    message ("-------------------------------------------------------------------------------------")
    message ("-- WILL COMPILE MPHF! (COMPILER VERSION IS HIGH ENOUGH) ")
    message ("-------------------------------------------------------------------------------------")
    if (SSE4_2_FOUND)
        set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS}    -DEMPHF_USE_POPCOUNT=1")
    endif()

    # for emphf
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        message ("-- EMPHF CLANG<3.8 FIX ")

        if (CMAKE_CXX_COMPILER_VERSION  VERSION_LESS 3.8)  # TODO perhaps might need to be revisited for Mac
            # my clang 3.7 (linux) doesn't like when i use stdlib=libc++, it refuses to compile because cannot found standard libs like <vector>
        else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        endif()

    endif ()

endif()

# In case we use an "old" version of clang with boost and c++0x, we have to skip rvalue usage
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	set (LIBRARY_COMPILE_DEFINITIONS  "${LIBRARY_COMPILE_DEFINITIONS} -DBOOST_NO_CXX11_RVALUE_REFERENCES=1")
endif()

# WARNING !!! For the moment, we need to remove some warnings (on Macos) due to use of offsetof macro on non Plain Old Data
set (LIBRARY_COMPILE_DEFINITIONS "${LIBRARY_COMPILE_DEFINITIONS} -Wno-invalid-offsetof") 

# We may have some custom compilation flags.
if (use_allocator)
    set (LIBRARY_COMPILE_DEFINITIONS "${LIBRARY_COMPILE_DEFINITIONS} -DGATB_USE_CUSTOM_ALLOCATOR") 
endif()

message("-- Options: ${LIBRARY_COMPILE_DEFINITIONS}")

################################################################################
#  DIRECTORIES MANAGEMENT 
################################################################################

set (LIBRARY_OUTPUT_PATH       ${PROJECT_BINARY_DIR}/lib/${CMAKE_BUILD_TYPE})
set (EXECUTABLE_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE})

################################################################################
#  EXPORT FOR THIRD PARTIES 
################################################################################

# In case you need to use external libraries, you can use the following variables:
#   EXTRALIBS      : names of the libraries to be used
#   EXTRALIBS_PATH : directory where the libraries are
#   EXTRALIBS_INC  : directory of the include files for the libraries 
# For instance, you can do the following:
#       cmake -DEXTRALIBS="somelib"  -DEXTRALIBS_PATH="/somewhere/libs"  -DEXTRALIBS_INC="/somewhere/inc" ..

IF (DEFINED EXTRALIBS)
    set (gatb-core-extra-libraries ${EXTRALIBS})
    set (gatb-core-extra-libraries ${gatb-core-extra-libraries} PARENT_SCOPE)
ENDIF()
IF (DEFINED EXTRALIBS_PATH)
    set (gatb-core-extra-libraries-path ${EXTRALIBS_PATH})
    set (gatb-core-extra-libraries-path ${gatb-core-extra-libraries-path} PARENT_SCOPE)
ENDIF()
IF (DEFINED EXTRALIBS_INC)
    set (gatb-core-extra-libraries-inc ${EXTRALIBS_INC})
    set (gatb-core-extra-libraries-inc ${gatb-core-extra-libraries-inc} PARENT_SCOPE)
ENDIF()

# We define the compilation flags used for compiling binary based on gatb core
set (gatb-core-flags ${LIBRARY_COMPILE_DEFINITIONS})

# We define the include directories used for linking binary based on gatb core
# Note that we need to add boost include dependency (not working otherwise with clang)
set (gatb-core-includes ${PROJECT_BINARY_DIR}/include  ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE} ${PROJECT_SOURCE_DIR}/src  ${PROJECT_SOURCE_DIR}/thirdparty ${gatb-core-extra-libraries-inc})

# We define the libraries used for linking binary based on gatb core
set (gatb-core-libraries   gatbcore-static  dl  pthread  z hdf5 ${gatb-core-extra-libraries})

# We define the directory where to find cmake helpers
set (gatb-core-cmake  ${CMAKE_MODULE_PATH})

# NOTE... we have to duplicate the variables for the other scopes (in particular for sub directories)
set (gatb-core-flags          ${gatb-core-flags}          PARENT_SCOPE)
set (gatb-core-includes       ${gatb-core-includes}       PARENT_SCOPE)
set (gatb-core-libraries      ${gatb-core-libraries}      PARENT_SCOPE)
set (gatb-core-cmake          ${gatb-core-cmake}          PARENT_SCOPE)

################################################################################
#  LIBRARY GENERATION 
################################################################################
ADD_SUBDIRECTORY(src)

# NOTE... we have to duplicate the variables for the other scopes (in particular for sub directories)
set (gatb-core-klist   ${gatb-core-klist}  PARENT_SCOPE)

################################################################################
#  UNIT TESTS GENERATION 
################################################################################
IF (DEFINED CPPUNIT_FOUND)
    IF (EXISTS "${PROJECT_SOURCE_DIR}/test")
        IF (NOT DEFINED GATB_CORE_EXCLUDE_TESTS)
            ADD_SUBDIRECTORY (test)
        ENDIF()
    ENDIF()
ENDIF()

################################################################################
#  TOOLS GENERATION 
################################################################################
IF (NOT DEFINED GATB_CORE_EXCLUDE_TOOLS)
    ADD_SUBDIRECTORY(tools)
ENDIF()

################################################################################
#  THIRD PARTY GENERATION (
################################################################################
ADD_SUBDIRECTORY(thirdparty)

################################################################################
#  DEPENDENCIES 
################################################################################
# we must be sure that hdf5 is built and installed before building gatb-core
ADD_DEPENDENCIES (gatbcore-static hdf5 hdf5_postbuild)

IF (DEFINED WITH_MPHF)
    ADD_DEPENDENCIES(gatbcore-static emphf_copyasis)
ENDIF()

################################################################################
#  DOCUMENTATION GENERATION 
################################################################################
IF (EXISTS "${PROJECT_SOURCE_DIR}/doc")
    ADD_SUBDIRECTORY(doc EXCLUDE_FROM_ALL)
ENDIF()

################################################################################
#  EXAMPLES GENERATION 
################################################################################
IF (EXISTS "${PROJECT_SOURCE_DIR}/examples")
    IF (NOT DEFINED GATB_CORE_EXCLUDE_EXAMPLES)
        ADD_SUBDIRECTORY(examples EXCLUDE_FROM_ALL)
    ENDIF()
ENDIF()

################################################################################
#  INSTALL 
################################################################################

IF (NOT DEFINED GATB_CORE_INSTALL_EXCLUDE)
    INSTALL (FILES ${PROJECT_SOURCE_DIR}/doc/misc/README.txt  DESTINATION . OPTIONAL)
    INSTALL (DIRECTORY ${PROJECT_SOURCE_DIR}/thirdparty/boost DESTINATION ./include)
ENDIF()

IF (NOT DEFINED GATB_CORE_EXCLUDE_EXAMPLES)
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
SET (CPACK_PACKAGE_FILE_NAME            "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-bin-${CMAKE_SYSTEM_NAME}")
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")
SET (CPACK_SOURCE_IGNORE_FILES          
    "^${PROJECT_SOURCE_DIR}/build/"  
    "^${PROJECT_SOURCE_DIR}/.project"
    "^${PROJECT_SOURCE_DIR}/.gitignore"
    "^${PROJECT_SOURCE_DIR}/doc/design"
    "^${PROJECT_SOURCE_DIR}/DELIVERY.md"
)



IF (NOT DEFINED GATB_CORE_INSTALL_EXCLUDE)
    # We include the module and get all the delivery targets
    include (Delivery)
ENDIF ()
