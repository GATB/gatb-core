################################################################################
# CPPUNIT
################################################################################

FIND_PATH (CPPUNIT_INCLUDE_DIR cppunit/extensions/HelperMacros.h
  $ENV{HOME}/.linuxbrew/include
  /local/include
  /usr/include
  /usr/local/Cellar/cppunit/1.12.1/include
  NO_DEFAULT_PATH
)

FIND_LIBRARY (CPPUNIT_LIBRARY cppunit
    ${CPPUNIT_INCLUDE_DIR}/../lib
)

# A little hack here... We check whether the library is reachable too.
if (NOT EXISTS "${CPPUNIT_INCLUDE_DIR}/../lib/libcppunit.a")
    message ("-- CppUnit: found headers (in ${CPPUNIT_INCLUDE_DIR}) but not the library...")
    SET (CPPUNIT_NO_LIB_FOUND 1)   
endif()

IF (CPPUNIT_INCLUDE_DIR)
    IF (CPPUNIT_LIBRARY)
        IF (NOT CPPUNIT_NO_LIB_FOUND)
            SET (CPPUNIT_FOUND "YES")
            SET (CPPUNIT_LIBRARIES ${CPPUNIT_LIBRARY})
            SET (CPPUNIT_LIBRARY_STATIC ${CPPUNIT_INCLUDE_DIR}/../lib/libcppunit.a)
        ENDIF()
    ENDIF (CPPUNIT_LIBRARY)
ENDIF (CPPUNIT_INCLUDE_DIR)

IF (DEFINED CPPUNIT_FOUND)
    message("-- CppUnit FOUND (${CPPUNIT_INCLUDE_DIR})")
ELSE()
    message("-- CppUnit NOT FOUND")
ENDIF()

