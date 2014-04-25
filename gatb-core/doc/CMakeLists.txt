
find_package (Doxygen)

# see http://www.cmake.org/pipermail/cmake/2008-May/021426.html

if (DOXYGEN_FOUND)

    CONFIGURE_FILE (${CMAKE_CURRENT_SOURCE_DIR}/doxygen/gatb.doxyfile  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile  @ONLY)
    
    ADD_CUSTOM_TARGET (doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    
    INSTALL (CODE "execute_process(COMMAND ${CMAKE_BUILD_TOOL} doc)")
    
    INSTALL (DIRECTORY ${CMAKE_BINARY_DIR}/doc/html DESTINATION doc)
    
endif (DOXYGEN_FOUND)


# HERE IS A COMMAND LINE TO COPY THE HTML DOCUMENTATION TO THE GFORGE SERVER  (from the build directory)
#    scp doc/html/* edrezen@scm.gforge.inria.fr:/home/groups/gatb-core/htdocs/