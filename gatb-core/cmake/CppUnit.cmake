################################################################################
# CPPUNIT
################################################################################

FIND_PATH (CPPUNIT_INCLUDE_DIR cppunit/extensions/HelperMacros.h
  /local/include
  /usr/include
  /usr/local/Cellar/cppunit/1.12.1/include
)

FIND_LIBRARY (CPPUNIT_LIBRARY cppunit
    ${CPPUNIT_INCLUDE_DIR}/../lib
)

IF (CPPUNIT_INCLUDE_DIR)
    IF (CPPUNIT_LIBRARY)
        SET (CPPUNIT_FOUND "YES")
        SET (CPPUNIT_LIBRARIES ${CPPUNIT_LIBRARY})
        SET (CPPUNIT_LIBRARY_STATIC ${CPPUNIT_INCLUDE_DIR}/../lib/libcppunit.a)
    ENDIF (CPPUNIT_LIBRARY)
ENDIF (CPPUNIT_INCLUDE_DIR)

IF (DEFINED CPPUNIT_FOUND)
    message("-- CppUnit FOUND")
ELSE()
    message("-- CppUnit NOT FOUND")
ENDIF()

