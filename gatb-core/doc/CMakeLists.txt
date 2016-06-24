
find_package (Doxygen)

# see http://www.cmake.org/pipermail/cmake/2008-May/021426.html

if (DOXYGEN_FOUND)

    CONFIGURE_FILE (${CMAKE_CURRENT_SOURCE_DIR}/doxygen/gatb.doxyfile  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile  @ONLY)
    
    ADD_CUSTOM_TARGET (doc 
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    )
    ADD_CUSTOM_TARGET (deploy-doc 
      COMMAND tar -zcf ${CMAKE_BINARY_DIR}/doc/doc.tgz -C ${CMAKE_BINARY_DIR}/doc/html .
      COMMAND scp ${CMAKE_BINARY_DIR}/doc/doc.tgz ${gatb-doc-user-login}@scm.gforge.inria.fr:/home/groups/gatb-core/
      COMMAND ssh ${gatb-doc-user-login}@scm.gforge.inria.fr "rm -rf /home/groups/gatb-core/htdocs/doc/api/*"
      COMMAND ssh ${gatb-doc-user-login}@scm.gforge.inria.fr "tar -xf /home/groups/gatb-core/doc.tgz --no-overwrite-dir -C /home/groups/gatb-core/htdocs/doc/api"
    )
    
    INSTALL (CODE "execute_process(COMMAND ${CMAKE_BUILD_TOOL} doc)")
    
    INSTALL (DIRECTORY ${CMAKE_BINARY_DIR}/doc/html DESTINATION doc)
    
endif (DOXYGEN_FOUND)


# HERE IS A COMMAND LINE TO COPY THE HTML DOCUMENTATION TO THE GFORGE SERVER  (from the build directory)
#    