################################################################################
# DELIVERY
################################################################################

# We get the 'package' and 'package_source' targets
include (CPack)

# We get the user name
SET (CPACK_USER_NAME  $ENV{USER})

# We set the server URI
SET (CPACK_SERVER_ADDRESS    "${CPACK_USER_NAME}@scm.gforge.inria.fr")
SET (CPACK_SERVER_DIR_BIN "/home/groups/${PROJECT_NAME}/versions/bin/")
SET (CPACK_SERVER_DIR_SRC "/home/groups/${PROJECT_NAME}/versions/src/")

# We define the name of the bin and src targets
SET (CPACK_URI_BIN "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME}.tar.gz")
SET (CPACK_URI_SRC "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-Source.tar.gz")

# We define the location where the bin and src targets have to be uploaded
SET (CPACK_UPLOAD_URI_BIN "${CPACK_SERVER_ADDRESS}:${CPACK_SERVER_DIR_BIN}")
SET (CPACK_UPLOAD_URI_SRC "${CPACK_SERVER_ADDRESS}:${CPACK_SERVER_DIR_SRC}")

# We add a custom target for delivery
add_custom_target (delivery 
    COMMAND echo "-----------------------------------------------------------"
    COMMAND echo "DELIVERY FOR ${PROJECT_NAME}, VERSION ${CPACK_PACKAGE_VERSION}"
    COMMAND echo "-----------------------------------------------------------"
    COMMAND make -j8 package package_source
    COMMAND scp ${CPACK_URI_BIN} ${CPACK_UPLOAD_URI_BIN}
    COMMAND scp ${CPACK_URI_SRC} ${CPACK_UPLOAD_URI_SRC}
)

# We add a custom target for delivery binaries
add_custom_target (delivery_bin 
    COMMAND echo "-----------------------------------------------------------"
    COMMAND echo "DELIVERY BIN FOR ${PROJECT_NAME}, VERSION ${CPACK_PACKAGE_VERSION}"
    COMMAND echo "-----------------------------------------------------------"
    COMMAND make -j8 package
    COMMAND scp ${CPACK_URI_BIN} ${CPACK_UPLOAD_URI_BIN}
)

# We add a custom target for delivery sources
add_custom_target (delivery_src
    COMMAND echo "-----------------------------------------------------------"
    COMMAND echo "DELIVERY SRC FOR ${PROJECT_NAME}, VERSION ${CPACK_PACKAGE_VERSION}"
    COMMAND echo "-----------------------------------------------------------"
    COMMAND make -j8 package_source
    COMMAND scp ${CPACK_URI_SRC} ${CPACK_UPLOAD_URI_SRC}
)

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

# We add a custom target for dumping existing deliveries
add_custom_target (delivery_dump 
    COMMAND echo ""
    COMMAND echo "-----------------------------------------------------------"
    COMMAND echo "EXISTING BINARIES..."
    COMMAND ssh  ${CPACK_SERVER_ADDRESS} "ls ${CPACK_SERVER_DIR_BIN}"
    COMMAND echo "-----------------------------------------------------------"
    COMMAND echo ""
    COMMAND echo "-----------------------------------------------------------"
    COMMAND echo "EXISTING SOURCES..."
    COMMAND ssh  ${CPACK_SERVER_ADDRESS} "ls ${CPACK_SERVER_DIR_SRC}"
    COMMAND echo "-----------------------------------------------------------"
    COMMAND echo ""
)
