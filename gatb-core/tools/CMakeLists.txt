################################################################################
#  EXAMPLES 
################################################################################

include_directories (${gatb-core-includes})

list (APPEND PROGRAMS dbgh5 dbgcheck bankgen)

FOREACH (program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${gatb-core-libraries})
ENDFOREACH (program)

################################################################################
#  INSTALLATION 
################################################################################
install (TARGETS dbgh5   DESTINATION bin)
