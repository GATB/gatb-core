################################################################################
#  EXAMPLES 
################################################################################

# We add the compilation options for the library
add_definitions (${gatb-core-flags})

include_directories (${gatb-core-includes})

list (APPEND PROGRAMS dbgh5 dbginfo dbgtopology bankgen)

FOREACH (program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${gatb-core-libraries})
ENDFOREACH (program)

################################################################################
#  INSTALLATION 
################################################################################
install (TARGETS dbgh5 dbginfo dbgtopology DESTINATION bin)
