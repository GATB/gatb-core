################################################################################
# CPPUNIT
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

IF (DEFINED CPPUNIT_FOUND)
    message("-- CppUnit FOUND")
ELSE()
    message("-- CppUnit NOT FOUND")
ENDIF()

