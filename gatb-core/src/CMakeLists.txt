################################################################################
#  CONFIGURATION FILE
################################################################################
configure_file (
    ${PROJECT_SOURCE_DIR}/src/gatb/system/api/config.hpp.in 
    ${PROJECT_BINARY_DIR}/include/gatb/system/api/config.hpp
)

install (FILES ${PROJECT_BINARY_DIR}/include/gatb/system/api/config.hpp DESTINATION include/gatb/system/api/)

################################################################################
#  HEADERS COPY
################################################################################
file (COPY ${PROJECT_SOURCE_DIR}/src/
    DESTINATION ${PROJECT_BINARY_DIR}/include
    FILES_MATCHING PATTERN "*.hpp"  PATTERN "*.tpp"  PATTERN "*.pri"
)

install (DIRECTORY ${PROJECT_SOURCE_DIR}/src/ DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.tpp" PATTERN "*.pri")

################################################################################
#  LIBRARY GENERATION 
################################################################################

file (GLOB_RECURSE  ComponentFiles  *)

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

# WARNING !!! For the moment, we need to remove some warnings (on Macos) due to use of offsetof macro on non Plain Old Data
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof") 

include_directories (${gatb-core-includes})

add_library    (gatbcore-static  STATIC  ${ComponentFiles} )
#add_library   (gatbcore-dynamic SHARED  ${ComponentFiles} )

set_target_properties  (gatbcore-static   PROPERTIES OUTPUT_NAME  gatbcore  clean_direct_output 1)
#set_target_properties (gatbcore-dynamic  PROPERTIES OUTPUT_NAME  gatbcore  clean_direct_output 1)

################################################################################
#  INSTALLATION 
################################################################################

# We install the libraries
install (TARGETS gatbcore-static DESTINATION lib)
#install (TARGETS gatbcore-dynamic DESTINATION lib)
