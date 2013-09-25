################################################################################
#  HDF5 GENERATION 
################################################################################

SET (HDF5_ENABLE_THREADSAFE       ON)

##########  MOMENTARY DEACTIVATED => CRASH ON MACOS TO BE INVESTIGATED  ##########
#SET (HDF5_BUILD_TOOLS             ON CACHE BOOL "Build HDF5 Tools")
#SET (CMAKE_EXE_LINKER_FLAGS       "-lpthread -lz")

SET (HDF5_EXTERNALLY_CONFIGURED     ON)
SET (HDF5_INSTALL_BIN_DIR           ${PROJECT_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE})
SET (HDF5_INSTALL_LIB_DIR           ${PROJECT_BINARY_DIR}/lib/${CMAKE_BUILD_TYPE})
SET (HDF5_INSTALL_INCLUDE_DIR       ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE})
SET (HDF5_INSTALL_DATA_DIR          ${PROJECT_BINARY_DIR}/share/${CMAKE_BUILD_TYPE})
SET (HDF5_INSTALL_CMAKE_DIR         ${PROJECT_BINARY_DIR}/share/${CMAKE_BUILD_TYPE})

# We don't want warnings from HDF5 compilation
set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS}  -w")
add_definitions (${COMPILE_DEFINITIONS})

# add HDF5 generation
ADD_SUBDIRECTORY (hdf5)
