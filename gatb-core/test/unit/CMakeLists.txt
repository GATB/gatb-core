################################################################################
#  CPPUNIT 
################################################################################
FIND_PATH (CPPUNIT_INCLUDE_DIR cppunit/extensions/HelperMacros.h
  /local/users/edrezen/local/include
  /usr/include
)

FIND_LIBRARY (CPPUNIT_LIBRARY cppunit
    ${CPPUNIT_INCLUDE_DIR}/../lib
)

IF (CPPUNIT_INCLUDE_DIR)
    IF (CPPUNIT_LIBRARY)
        SET (CPPUNIT_FOUND "YES")
        SET (CPPUNIT_LIBRARIES ${CPPUNIT_LIBRARY})
    ENDIF (CPPUNIT_LIBRARY)
ENDIF (CPPUNIT_INCLUDE_DIR)

################################################################################
#  COMPONENT 
################################################################################

set (LIBRARY_OUTPUT_PATH    lib/${CMAKE_BUILD_TYPE})

file (GLOB_RECURSE  GatbToolsTestFiles  src/*)

include_directories (${org.inria.gatb.tools_SOURCE_DIR}/src  ${org.inria.gatb.tools_SOURCE_DIR}/test/unit/src  ${CPPUNIT_INCLUDE_DIR})

add_executable        (GatbToolsTest    ${GatbToolsTestFiles})
target_link_libraries (GatbToolsTest     GatbTools-static  ${CPPUNIT_LIBRARY}  dl pthread z)

################################################################################
#  LAUNCH 
################################################################################

# We add an optional target for launching unit tests
add_custom_target (test GatbToolsTest ${org.inria.gatb.tools_SOURCE_DIR}/test/db  unittests.xml)
