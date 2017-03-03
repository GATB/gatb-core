################################################################################
#  CPPUNIT 
################################################################################

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

file (GLOB_RECURSE  GatbCoreTestFiles  src/*)

include_directories (${gatb-core-includes}  ${gatb-core_SOURCE_DIR}/test/unit/src  ${CPPUNIT_INCLUDE_DIR})

# We add the path for extra libraries
link_directories (${gatb-core-extra-libraries-path})

add_executable        (gatb-core-cppunit    ${GatbCoreTestFiles})
target_link_libraries (gatb-core-cppunit    ${gatb-core-libraries} ${CPPUNIT_LIBRARIES})

#################################################
#  LAUNCH 
#################################################

# We add an optional target for launching unit tests
add_custom_target (gatb_test gatb-core-cppunit  "All Tests"  ${gatb-core_SOURCE_DIR}/test/db  unittests.xml)

#################################################
#  INSTALL 
#################################################
install (TARGETS gatb-core-cppunit   DESTINATION test)
install (DIRECTORY ${gatb-core_SOURCE_DIR}/test/db DESTINATION test/)
    
