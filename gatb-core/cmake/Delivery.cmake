################################################################################
# DELIVERY
################################################################################

# We get the 'package' and 'package_source' targets
include (CPack)

# We get the user name
IF (NOT CPACK_USER_NAME)
    SET (CPACK_USER_NAME  $ENV{USER})
ENDIF (NOT CPACK_USER_NAME)

# We get the date
GetCurrentDateShort (CPACK_DATE) 

# We set the name of the versions file.
SET (CPACK_VERSIONS_FILENAME  "versions.txt")

# We may have to set (if not defined) the CPACK_GFORGE_PROJECT_NAME
IF (NOT CPACK_GFORGE_PROJECT_NAME)
    SET (CPACK_GFORGE_PROJECT_NAME  ${PROJECT_NAME})
ENDIF (NOT CPACK_GFORGE_PROJECT_NAME)

# We set the server URI
SET (CPACK_SERVER_ADDRESS   "${CPACK_USER_NAME}@scm.gforge.inria.fr")
SET (CPACK_SERVER_DIR       "/home/groups/${CPACK_GFORGE_PROJECT_NAME}/htdocs/versions/")
SET (CPACK_SERVER_DIR_BIN   "${CPACK_SERVER_DIR}/bin/")
SET (CPACK_SERVER_DIR_SRC   "${CPACK_SERVER_DIR}/src/")
SET (CPACK_SERVER_VERSIONS  "${CPACK_SERVER_DIR}/${CPACK_VERSIONS_FILENAME}")

# We define the name of the bin and src targets
SET (CPACK_URI_BIN "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME}.tar.gz")
SET (CPACK_URI_SRC "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-Source.tar.gz")

# We define the location where the bin and src targets have to be uploaded
SET (CPACK_UPLOAD_URI_BIN  "${CPACK_SERVER_ADDRESS}:${CPACK_SERVER_DIR_BIN}")
SET (CPACK_UPLOAD_URI_SRC  "${CPACK_SERVER_ADDRESS}:${CPACK_SERVER_DIR_SRC}")
SET (CPACK_UPLOAD_VERSIONS "${CPACK_SERVER_ADDRESS}:${CPACK_SERVER_VERSIONS}")

# We set the text holding all the information about the delivery.
SET (CPACK_INFO_BIN ${CPACK_GFORGE_PROJECT_NAME} bin ${CMAKE_PROJECT_NAME} ${CPACK_PACKAGE_VERSION} ${CPACK_DATE} ${CPACK_SYSTEM_NAME} ${CPACK_USER_NAME} ${CPACK_URI_BIN})
SET (CPACK_INFO_SRC ${CPACK_GFORGE_PROJECT_NAME} src ${CMAKE_PROJECT_NAME} ${CPACK_PACKAGE_VERSION} ${CPACK_DATE} all                  ${CPACK_USER_NAME} ${CPACK_URI_SRC})


################################################################################
# MAIN TARGET 
################################################################################

# We add a custom target for delivery
add_custom_target (delivery     

    DEPENDS delivery_bin  delivery_src 
    
    COMMAND echo "-----------------------------------------------------------"
    COMMAND echo "DELIVERY FOR ${CPACK_GFORGE_PROJECT_NAME}, VERSION ${CPACK_PACKAGE_VERSION}"
    COMMAND echo "-----------------------------------------------------------"

    # We dump the known versions
    COMMAND make delivery_dump
)

################################################################################
# TARGETS 'bin'
################################################################################

# We add a custom target for delivery binaries
add_custom_target (delivery_bin 

    # We get the versions.txt file from the server
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/delivery.sh  "BIN" ${CPACK_GFORGE_PROJECT_NAME} ${CPACK_PACKAGE_VERSION} ${CPACK_UPLOAD_VERSIONS} ${CPACK_VERSIONS_FILENAME}  \"${CPACK_INFO_BIN}\"  ${CPACK_URI_BIN}   ${CPACK_UPLOAD_URI_BIN}
)

################################################################################
# TARGETS 'src'
################################################################################

# We add a custom target for delivery sources
add_custom_target (delivery_src 

    # We get the versions.txt file from the server
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/delivery.sh  "SRC" ${CPACK_GFORGE_PROJECT_NAME} ${CPACK_PACKAGE_VERSION} ${CPACK_UPLOAD_VERSIONS} ${CPACK_VERSIONS_FILENAME}  \"${CPACK_INFO_SRC}\"  ${CPACK_URI_SRC}   ${CPACK_UPLOAD_URI_SRC}
)

################################################################################
# TARGET 'help'
################################################################################

# We add a custom target for delivery sources
add_custom_target (delivery_help
    COMMAND echo "-----------------------------------------------------------"
    COMMAND echo "DELIVERY TARGETS"
    COMMAND echo "-----------------------------------------------------------"
    COMMAND echo "delivery:      build a targz for binaries and sources and upload them to GForge"
    COMMAND echo "delivery_bin:  build a targz for binaries and upload it to GForge"
    COMMAND echo "delivery_src:  build a targz for sources and upload it to GForge"
    COMMAND echo "delivery_dump: dump existing releases on GForge"
    COMMAND echo ""
)

################################################################################
# TARGET 'dump'
################################################################################

# We add a custom target for dumping existing deliveries
add_custom_target (delivery_dump

    # We get the versions.txt file from the server
    COMMAND scp ${CPACK_UPLOAD_VERSIONS} ${CPACK_VERSIONS_FILENAME}

    # We dump the versions file.
    COMMAND echo ""
    COMMAND echo "-------------------------------------------------------------------------------------------------"
    COMMAND echo "LIST OF DELIVERIES FOR " ${CPACK_GFORGE_PROJECT_NAME}
    COMMAND echo "-------------------------------------------------------------------------------------------------"
    COMMAND cat ${CPACK_VERSIONS_FILENAME}
    COMMAND echo "-------------------------------------------------------------------------------------------------"
    COMMAND echo ""
)
