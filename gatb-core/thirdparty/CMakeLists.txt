################################################################################
#  HDF5 GENERATION 
################################################################################

#SET (HDF5_ENABLE_THREADSAFE       ON)
#SET (H5_HAVE_THREADSAFE 1)

##########  MOMENTARY DEACTIVATED => CRASH ON MACOS TO BE INVESTIGATED  ##########
SET (HDF5_BUILD_TOOLS             ON CACHE BOOL "Build HDF5 Tools")
#SET (CMAKE_EXE_LINKER_FLAGS       "-lpthread -lz")

SET (HDF5_EXTERNALLY_CONFIGURED     ON)

#SET (HDF5_INSTALL_BIN_DIR           ${PROJECT_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE})
#SET (HDF5_INSTALL_LIB_DIR           ${PROJECT_BINARY_DIR}/lib/${CMAKE_BUILD_TYPE})
SET (HDF5_INSTALL_BIN_DIR           bin)
SET (HDF5_INSTALL_LIB_DIR           lib)

SET (HDF5_INSTALL_INCLUDE_DIR       ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/hdf5)
SET (HDF5_INSTALL_DATA_DIR          ${PROJECT_BINARY_DIR}/share/${CMAKE_BUILD_TYPE})
SET (HDF5_INSTALL_CMAKE_DIR         ${PROJECT_BINARY_DIR}/share/${CMAKE_BUILD_TYPE})

IF (NOT DEFINED GATB_CORE_INSTALL_EXCLUDE)
    SET (HDF5_EXPORTED_TARGETS          "gatb-hdf5")
ENDIF()

IF (NOT DEFINED GATB_CORE_EXCLUDE_HDF5_ZLIB)
	OPTION (HDF5_ENABLE_Z_LIB_SUPPORT "Enable Zlib Filters" ON)
ENDIF()

# We don't want warnings from HDF5 compilation
set (COMPILE_DEFINITIONS  "${COMPILE_DEFINITIONS}  -w")
add_definitions (${COMPILE_DEFINITIONS})

# add HDF5 generation
ADD_SUBDIRECTORY (hdf5)

# We add a custom target for copying header files.
add_custom_target (hdf5_postbuild ALL)

# We build the output directory
add_custom_command (TARGET hdf5_postbuild POST_BUILD  COMMAND ${CMAKE_COMMAND} -E make_directory ${HDF5_INSTALL_INCLUDE_DIR})

# We define all the header files to be copied
file (GLOB headerfiles ${PROJECT_SOURCE_DIR}/thirdparty/hdf5/src/*.h ${PROJECT_BINARY_DIR}/thirdparty/hdf5/H5pubconf.h)

# We copy each header file
foreach (header ${headerfiles})
  add_custom_command (TARGET hdf5_postbuild POST_BUILD  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${header} ${HDF5_INSTALL_INCLUDE_DIR}  )
endforeach()

################################################################################
#  EMPHF GENERATION 
################################################################################

# TODO: I actually just mimicked hdf5, but I'm not so sure anymore whether
# this is useful. MPHF is just #include'd, perhaps we can just not copy it? in
# that case just delete this code

IF (DEFINED use_mphf)

    ADD_SUBDIRECTORY (emphf)

    # We add a custom target for copying emphf as it is.
    add_custom_target (emphf_copyasis ALL)

    SET (EMPHF_INSTALL_INCLUDE_DIR
        ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/emphf)

    # We define all the header files to be copied
    file (GLOB headerfiles ${PROJECT_SOURCE_DIR}/thirdparty/emphf/*.hpp
        ${PROJECT_SOURCE_DIR}/thirdparty/emphf/*.cpp)

    # create emphf directory
    add_custom_command (TARGET emphf_copyasis COMMAND ${CMAKE_COMMAND} -E
        make_directory ${EMPHF_INSTALL_INCLUDE_DIR})

    # We copy each header file
    foreach (header ${headerfiles})
      add_custom_command (TARGET emphf_copyasis 
          COMMAND ${CMAKE_COMMAND} -E copy_if_different ${header} ${EMPHF_INSTALL_INCLUDE_DIR}  
        )
    endforeach()
ENDIF() # WITH_MPHF

################################################################################
#  INSTALL 
################################################################################
IF (NOT DEFINED GATB_CORE_INSTALL_EXCLUDE)
    INSTALL (DIRECTORY ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/hdf5 DESTINATION include)
    IF (DEFINED WITH_MPHF)
        INSTALL (DIRECTORY ${PROJECT_BINARY_DIR}/include/${CMAKE_BUILD_TYPE}/emphf DESTINATION include)
    ENDIF()
ENDIF ()
