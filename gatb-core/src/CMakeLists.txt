################################################################################
#  KMER SIZE DEFINITION 
################################################################################

# We define several sizes of kmers
SET (KSIZE_1  "32")
SET (KSIZE_2  "64")
SET (KSIZE_3  "96")
SET (KSIZE_4 "128")

# Users can overwrite the default fourth value
IF (kmax)
    SET (KSIZE_4 ${kmax})
ENDIF()

MESSAGE ("-- SUPPORTED KMER SIZES ARE ${KSIZE_1} ${KSIZE_2} ${KSIZE_3} ${KSIZE_4}")

################################################################################
#  CONFIGURATION FILE
################################################################################
configure_file (
    ${PROJECT_SOURCE_DIR}/src/gatb/system/api/config.hpp.in 
    ${PROJECT_BINARY_DIR}/include/gatb/system/api/config.hpp
)

################################################################################
#  LIBRARY GENERATION 
################################################################################

file (GLOB_RECURSE  ComponentFiles  *)

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

# WARNING !!! For the moment, we need to remove some warnings (on Macos) due to use of offsetof macro on non Plain Old Data
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof -Wno-ambiguous-member-template") 

include_directories (${gatb-core-includes})

add_library    (gatbcore-static  STATIC  ${ComponentFiles} )
#add_library   (gatbcore-dynamic SHARED  ${ComponentFiles} )

set_target_properties  (gatbcore-static   PROPERTIES OUTPUT_NAME  gatbcore  clean_direct_output 1)
#set_target_properties (gatbcore-dynamic  PROPERTIES OUTPUT_NAME  gatbcore  clean_direct_output 1)

################################################################################
#  INSTALLATION 
################################################################################

# We install the libraries
IF (NOT DEFINED GATB_CORE_INSTALL_EXCLUDE)
    install (TARGETS gatbcore-static DESTINATION lib)
    install (DIRECTORY ${PROJECT_SOURCE_DIR}/src/ DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.tpp" PATTERN "*.pri")
    install (FILES ${PROJECT_BINARY_DIR}/include/gatb/system/api/config.hpp DESTINATION include/gatb/system/api/)
ENDIF ()
