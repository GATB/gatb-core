cmake_minimum_required (VERSION 3.1.0)
PROJECT (HDF5_TOOLS_H5DUMP)

#-----------------------------------------------------------------------------
# Setup include Directories
#-----------------------------------------------------------------------------
INCLUDE_DIRECTORIES (${HDF5_TOOLS_SRC_DIR}/lib)

# --------------------------------------------------------------------
# Add the h5dump executables
# --------------------------------------------------------------------
add_executable (gatb-h5dump
    ${HDF5_TOOLS_H5DUMP_SOURCE_DIR}/h5dump.c
    ${HDF5_TOOLS_H5DUMP_SOURCE_DIR}/h5dump_ddl.c
    ${HDF5_TOOLS_H5DUMP_SOURCE_DIR}/h5dump_xml.c
)
TARGET_NAMING (gatb-h5dump STATIC)
TARGET_C_PROPERTIES (gatb-h5dump STATIC " " " ")
target_link_libraries (gatb-h5dump  ${HDF5_TOOLS_LIB_TARGET} ${HDF5_LIB_TARGET})
set_target_properties (gatb-h5dump PROPERTIES FOLDER tools)
set_global_variable (HDF5_UTILS_TO_EXPORT "${HDF5_UTILS_TO_EXPORT};gatb-h5dump")

# (rayan) attempt at static linking, not sure if it works properly
if (static)
    set_target_properties(gatb-h5dump PROPERTIES LINK_SEARCH_START_STATIC 1)
    set_target_properties(gatb-h5dump PROPERTIES LINK_SEARCH_END_STATIC 1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    if (APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++") #removed -static because doesn't work on mac's
    else()
        set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++ -static")
    endif()
   message ("STATIC BINARIES for gatb-h5dump")
else()
   message ("DYNAMIC BINARIES for gatb-h5dump")
endif()


set (H5_DEP_EXECUTABLES gatb-h5dump)

#if (BUILD_TESTING)
#  # --------------------------------------------------------------------
#  # Add the h5dump test executable
#  # --------------------------------------------------------------------
#  if (HDF5_BUILD_GENERATORS)
#    add_executable (h5dumpgentest ${HDF5_TOOLS_H5DUMP_SOURCE_DIR}/h5dumpgentest.c)
#    TARGET_NAMING (h5dumpgentest STATIC)
#    TARGET_C_PROPERTIES (h5dumpgentest STATIC " " " ")
#    target_link_libraries (h5dumpgentest ${HDF5_LIB_TARGET} ${HDF5_TOOLS_LIB_TARGET})
#    set_target_properties (h5dumpgentest PROPERTIES FOLDER generator/tools)
#
#   #add_test (NAME h5dumpgentest COMMAND $<TARGET_FILE:h5dumpgentest>)
#  endif (HDF5_BUILD_GENERATORS)
#
#  include (CMakeTests.cmake)
#
#  include (CMakeTestsPBITS.cmake)
#
#  include (CMakeTestsXML.cmake)
#endif (BUILD_TESTING)

##############################################################################
##############################################################################
###           I N S T A L L A T I O N                                      ###
##############################################################################
##############################################################################

#-----------------------------------------------------------------------------
# Rules for Installation of tools using make Install target
#-----------------------------------------------------------------------------

#INSTALL_PROGRAM_PDB (gatb-h5dump ${HDF5_INSTALL_BIN_DIR} toolsapplications)

install (
    TARGETS
        gatb-h5dump
    EXPORT
        ${HDF5_EXPORTED_TARGETS}
    RUNTIME DESTINATION ${HDF5_INSTALL_BIN_DIR} COMPONENT toolsapplications
)
