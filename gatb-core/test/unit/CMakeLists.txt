################################################################################
#  CPPUNIT 
################################################################################

file (GLOB_RECURSE  GatbToolsTestFiles  src/*)

include_directories (${gatb-core-includes}  ${org.inria.gatb.tools_SOURCE_DIR}/test/unit/src  ${CPPUNIT_INCLUDE_DIR})

add_executable        (GatbToolsTest    ${GatbToolsTestFiles})
target_link_libraries (GatbToolsTest    ${gatb-core-libraries} ${CPPUNIT_LIBRARY})

################################################################################
#  LAUNCH 
################################################################################

# We add an optional target for launching unit tests
add_custom_target (test GatbToolsTest  "All Tests"  ${org.inria.gatb.tools_SOURCE_DIR}/test/db  unittests.xml)
