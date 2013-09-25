################################################################################
#  EXAMPLES 
################################################################################

include_directories (${gatb-core-includes})

list (APPEND PROGRAMS dbgh5 dbgcheck)

FOREACH (program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${gatb-core-libraries})
ENDFOREACH (program)
