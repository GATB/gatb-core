################################################################################
#  EXAMPLES 
################################################################################

# We add the compilation options for the library
add_definitions (${gatb-core-flags} " -Wno-ambiguous-member-template")

include_directories (${gatb-core-includes})

list (APPEND PROGRAMS dbgh5 dbginfo)

FOREACH (program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${gatb-core-libraries})
ENDFOREACH (program)

################################################################################
#  INSTALLATION 
################################################################################
IF (NOT DEFINED GATB_CORE_INSTALL_EXCLUDE)
    install (TARGETS dbgh5 dbginfo DESTINATION bin)
ENDIF ()
