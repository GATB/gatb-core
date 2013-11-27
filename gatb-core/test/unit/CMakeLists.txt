################################################################################
#  CPPUNIT 
################################################################################

file (GLOB_RECURSE  GatbCoreTestFiles  src/*)

include_directories (${gatb-core-includes}  ${org.inria.gatb.tools_SOURCE_DIR}/test/unit/src  ${CPPUNIT_INCLUDE_DIR})

add_executable        (GatbCoreTest    ${GatbCoreTestFiles})
target_link_libraries (GatbCoreTest    ${gatb-core-libraries} ${CPPUNIT_LIBRARY})

################################################################################
#  LAUNCH 
################################################################################

# We add an optional target for launching unit tests
add_custom_target (test GatbCoreTest  "All Tests"  ${org.inria.gatb.tools_SOURCE_DIR}/test/db  unittests.xml)
