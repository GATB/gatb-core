################################################################################
#  KMER SIZE DEFINITION 
################################################################################

# We define several sizes of kmers
SET (KSIZE_1  "32")
SET (KSIZE_2  "64")
SET (KSIZE_3  "96")
SET (KSIZE_4 "128")
SET (KSIZE_5 "160")
SET (KSIZE_6 "192")
SET (KSIZE_7 "224")
SET (KSIZE_8 "256")

# Users can overwrite these values if needed
IF (k1)
    SET (KSIZE_1 ${k1})
ENDIF()

IF (k2)
    SET (KSIZE_2 ${k2})
ENDIF()

IF (k3)
    SET (KSIZE_3 ${k3})
ENDIF()

IF (k4)
    SET (KSIZE_4 ${k4})
ENDIF()

IF (k5)
    SET (KSIZE_5 ${k5})
ENDIF()

IF (k6)
    SET (KSIZE_6 ${k6})
ENDIF()

IF (k7)
    SET (KSIZE_7 ${k7})
ENDIF()

IF (k8)
    SET (KSIZE_8 ${k8})
ENDIF()



MESSAGE ("-- OPTIMIZED KMER SIZES INTERVALS ARE ${KSIZE_1} ${KSIZE_2} ${KSIZE_3} ${KSIZE_4} ${KSIZE_5} ${KSIZE_6} ${KSIZE_7} ${KSIZE_8} <-- max supported kmer size without recompilation")

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
